// module2/js/gc-balancer.js

import { CONFIG } from './config.js';

export class GCBalancer {
    constructor(codonTable) {
        this.codonTable = codonTable;
    }

    balance(dna, targetGC = CONFIG.gc.defaultTarget, windowSize = CONFIG.gc.defaultWindowSize, tolerance = CONFIG.gc.tolerance) {
        let currentDna = dna;
        const codons = [];
        for (let i = 0; i < currentDna.length; i += 3) {
            codons.push(currentDna.substring(i, i + 3));
        }

        let iterations = 0;
        const maxIterations = CONFIG.gc.maxIterations;

        while (iterations < maxIterations) {
            let changed = false;
            for (let i = 0; i <= currentDna.length - windowSize; i += 3) {
                const window = currentDna.substring(i, i + windowSize);
                const currentGC = this.calculateGC(window);

                if (Math.abs(currentGC - targetGC) > tolerance) {
                    const startCodonIdx = Math.floor(i / 3);
                    const endCodonIdx = Math.floor((i + windowSize) / 3);

                    for (let j = startCodonIdx; j < endCodonIdx; j++) {
                        const aa = this.codonTable.getAAForCodon(codons[j]);
                        // Filter for synonymous codons that are high-quality (freq > 50% of max for that AA)
                        const allSynonymous = this.codonTable.getCodonsForAA(aa);
                        const maxFreq = allSynonymous[0].freq;
                        const minFreq = maxFreq * CONFIG.gc.minFreqPercentile;
                        const synonymous = allSynonymous.filter(c => c.freq >= minFreq);

                        if (synonymous.length > 1) {
                            let bestCodon = codons[j];
                            let bestGC = currentGC;

                            for (const s of synonymous) {
                                if (s.codon === codons[j]) continue;
                                
                                const newCodons = [...codons];
                                newCodons[j] = s.codon;
                                const newWindow = newCodons.slice(startCodonIdx, endCodonIdx).join('');
                                const newGC = this.calculateGC(newWindow);

                                if (Math.abs(newGC - targetGC) < Math.abs(bestGC - targetGC)) {
                                    bestGC = newGC;
                                    bestCodon = s.codon;
                                }
                            }

                            if (bestCodon !== codons[j]) {
                                codons[j] = bestCodon;
                                changed = true;
                            }
                        }
                    }
                }
            }

            currentDna = codons.join('');
            if (!changed) break;
            iterations++;
        }

        return currentDna;
    }

    calculateGC(sequence) {
        if (sequence.length === 0) return 0;
        const gcCount = (sequence.match(/[GC]/gi) || []).length;
        return gcCount / sequence.length;
    }
}
