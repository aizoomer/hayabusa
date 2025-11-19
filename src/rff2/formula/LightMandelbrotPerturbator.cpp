//
// Created by Merutilm on 2025-05-10.
// Created by Super Fractal on 2025-11-19.
//

#include "LightMandelbrotPerturbator.h"

#include <cmath>
#include "Perturbator.h"

namespace merutilm::rff2 {

    LightMandelbrotPerturbator::LightMandelbrotPerturbator(ParallelRenderState &state, const FractalAttribute &calc,
                                                           const double dcMax, const int exp10,
                                                           const uint64_t initialPeriod,
                                                           ApproxTableCache &tableRef,
                                                           std::function<void(uint64_t)> &&actionPerRefCalcIteration,
                                                           std::function<void(uint64_t, double)> &&
                                                           actionPerCreatingTableIteration,
                                                           const bool arbitraryPrecisionFPGBn,
                                                           std::unique_ptr<LightMandelbrotReference> reusedReference,
                                                           std::unique_ptr<LightMPATable> reusedTable,
                                                           const double offR,
                                                           const double offI) : MandelbrotPerturbator(state, calc),
                                                                                dcMax(dcMax), offR(offR), offI(offI) {
        if (reusedReference == nullptr) {
            reference = LightMandelbrotReference::createReference(state, calc, exp10, initialPeriod, dcMax,
                                                                  arbitraryPrecisionFPGBn,
                                                                  std::move(actionPerRefCalcIteration));
        } else {
            reference = std::move(reusedReference);
        }

        if (reference == Constants::NullPointer::PROCESS_TERMINATED_REFERENCE) {
            return;
        }

        if (reusedTable == nullptr) {
            table = std::make_unique<LightMPATable>(state, *reference, &calc.mpaAttribute, dcMax,
                                                    tableRef,
                                                    std::move(actionPerCreatingTableIteration));
        } else {
            table = std::move(reusedTable);
        }
    }

    double LightMandelbrotPerturbator::iterate(const dex &dcr, const dex &dci) const {
        if (state.interruptRequested()) return 0.0;

        // 定数をローカル変数へ
        const double dcr1 = static_cast<double>(dcr) + offR;
        const double dci1 = static_cast<double>(dci) + offI;

        uint64_t iteration = 0;
        uint64_t refIteration = 0;
        int absIteration = 0;
        const uint64_t maxRefIteration = reference->longestPeriod();

        double dzr = 0; // delta z
        double dzi = 0;
        double zr = 0; // z
        double zi = 0;

        double cd = 0;
        double pd = 0;
        
        const bool isAbs = calc.absoluteIterationMode;
        const uint64_t maxIteration = calc.maxIteration;
        const float bailout = calc.bailout;
        const float bailout2 = bailout * bailout;
        const int exitCheckInterval = Constants::Fractal::EXIT_CHECK_INTERVAL;

        // --- Optimization: Access raw pointers directly to avoid indirection overhead ---
        // std::vector等のデータへの直接ポインタを取得（実装依存ですが、通常operator[]より高速）
        // referenceクラスの実装が見えないため、安全にoperator[]を使う形を維持しつつ、
        // compressorとrefオブジェクトへのポインタをキャッシュします。
        const auto* refObj = reference.get();
        const auto* mpaTable = table.get();
        
        // 最初の参照軌道をロード
        uint64_t index = ArrayCompressor::compress(refObj->compressor, refIteration);
        // Note: refReal/refImagがpublicメンバ変数であると仮定（元のコードに基づく）
        // 可能なら refReal.data() を取得したいが、std::vectorか不明なためそのまま使用
        // ループ内での間接参照を減らすため、現在値をキャッシュ
        double curRefR = refObj->refReal[index];
        double curRefI = refObj->refImag[index];

        // 中断チェック用カウンタ（剰余演算の除去）
        int checkCounter = exitCheckInterval;

        while (iteration < maxIteration) {
            // MPA Optimization
            // mpaTableがnullでない場合のみチェック。分岐予測によりコストは低い。
            if (mpaTable) {
                if (const LightPA *mpaPtr = mpaTable->lookup(refIteration, dzr, dzi); mpaPtr != nullptr) {
                    const LightPA &mpa = *mpaPtr;
                    // MPAによるスキップ計算
                    const double dzr1 = mpa.anr * dzr - mpa.ani * dzi + mpa.bnr * dcr1 - mpa.bni * dci1;
                    const double dzi1 = mpa.anr * dzi + mpa.ani * dzr + mpa.bnr * dci1 + mpa.bni * dcr1;

                    dzr = dzr1;
                    dzi = dzi1;

                    const uint32_t skip = mpa.skip;
                    iteration += skip;
                    refIteration += skip; // ここでrefIterationが大きく進む可能性がある
                    absIteration++;       // 元コードではskip時もabsIterationは+1のみ

                    if (iteration >= maxIteration) {
                        return static_cast<double>(isAbs ? absIteration : maxIteration);
                    }
                    
                    // MPAスキップ後、参照軌道のキャッシュを更新する必要がある
                    index = ArrayCompressor::compress(refObj->compressor, refIteration);
                    curRefR = refObj->refReal[index];
                    curRefI = refObj->refImag[index];
                    
                    // ここで continue するとループ条件チェックへ戻る
                    continue;
                }
            }

            // --- Perturbation Logic ---
            
            // 元コード: if (refIteration != maxRefIteration) ブロック内
            // maxRefIterationに到達している場合は摂動計算をスキップする仕様と見受けられますが、
            // 通常のマンデルブロセットでは稀なケースです。
            if (refIteration != maxRefIteration) {
                // Optimization: 2.0倍は加算の方が速い場合がある (コンパイラ最適化に委ねても良いが明示)
                // z = 2*Z*z + z^2 + c
                // Real: 2(Zr*zr - Zi*zi) + (zr^2 - zi^2) + cr -> (2Zr + zr)*zr - (2Zi + zi)*zi + cr
                const double twoZr_p_dzr = curRefR + curRefR + dzr;
                const double twoZi_p_dzi = curRefI + curRefI + dzi;
                
                const double zr2 = twoZr_p_dzr * dzr - twoZi_p_dzi * dzi + dcr1;
                const double zi2 = twoZr_p_dzr * dzi + twoZi_p_dzi * dzr + dci1;

                dzr = zr2;
                dzi = zi2;

                ++refIteration;
                ++iteration;
                ++absIteration;
                
                // 次の参照軌道を取得 (ArrayCompressor呼び出しをここに移動)
                // indexはループスコープ外の変数を利用
                index = ArrayCompressor::compress(refObj->compressor, refIteration);
                curRefR = refObj->refReal[index];
                curRefI = refObj->refImag[index];
            }
            
            // 現在の z = Ref + delta
            zr = curRefR + dzr;
            zi = curRefI + dzi;

            // Glitch / Validity check
            if (zi == 0.0 && zr < 0.25 && zr >= -2.0) {
                return static_cast<double>(maxIteration);
            }

            pd = cd;
            cd = zr * zr + zi * zi;

            // Rebasing / Reference wrapping logic
            if (refIteration == maxRefIteration || cd < dzr * dzr + dzi * dzi) {
                refIteration = 0;
                dzr = zr;
                dzi = zi;
                
                // 参照軌道をリセットしたため、キャッシュも更新
                index = ArrayCompressor::compress(refObj->compressor, 0);
                curRefR = refObj->refReal[index];
                curRefI = refObj->refImag[index];
            }

            if (cd > bailout2) {
                break;
            }

            // Optimization: Modulo removal
            if (--checkCounter == 0) {
                if (state.interruptRequested()) return 0.0;
                checkCounter = exitCheckInterval;
            }
        }

        if (isAbs) {
            return static_cast<double>(absIteration);
        }

        if (iteration >= maxIteration) {
            return static_cast<double>(maxIteration);
        }

        pd = std::sqrt(pd);
        cd = std::sqrt(cd);

        return getDoubleValueIteration(iteration, pd, cd, calc.decimalizeIterationMethod, bailout);
    }


    std::unique_ptr<LightMandelbrotPerturbator> LightMandelbrotPerturbator::reuse(
        const FractalAttribute &calc, const double dcMax, ApproxTableCache &tableRef) {

        const int exp10 = logZoomToExp10(calc.logZoom);
        double offR = 0;
        double offI = 0;
        uint64_t longestPeriod = 1;
        std::unique_ptr<LightMandelbrotReference> reusedReference = nullptr;

        if (reference == Constants::NullPointer::PROCESS_TERMINATED_REFERENCE) {
            MessageBox(nullptr, "Please do not try to use PROCESS-TERMINATED Reference.", "Warning",
                       MB_OK | MB_ICONWARNING);
        } else {
            fp_complex_calculator centerOffset = calc.center.edit(exp10);
            centerOffset -= reference->center.edit(exp10);
            offR = centerOffset.getReal().double_value();
            offI = centerOffset.getImag().double_value();
            longestPeriod = reference->longestPeriod();
            reusedReference = std::move(reference);
        }

        return std::make_unique<LightMandelbrotPerturbator>(state, calc, dcMax, exp10, longestPeriod, tableRef,
                                                            [](uint64_t) {}, [](uint64_t, double) {}, 
                                                            false, std::move(reusedReference),
                                                            std::move(table), offR, offI);
    }
}
