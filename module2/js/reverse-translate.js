// module2/js/reverse-translate.js

export class ReverseTranslator {
    constructor(codonTable) {
        this.codonTable = codonTable;
    }

    translate(protein, strategy, targetGC = 0.52) {
        switch (strategy) {
            case 'most_frequent':
                return this.mostFrequent(protein, targetGC);
            case 'weighted_random':
                return this.weightedRandom(protein, targetGC);
            case 'balanced':
                return this.balanced(protein, targetGC);
            default:
                throw new Error(`Unknown strategy: ${strategy}`);
        }
    }

    getCodonGC(codon) {
        return (codon.match(/[GC]/gi) || []).length / 3;
    }

    getGCWeight(codon, currentGC, targetGC) {
        const codonGC = this.getCodonGC(codon);
        // If we are below target, favor high GC codons. If above, favor low GC.
        const diff = targetGC - currentGC;
        // Exponential weight: the further we are from target, the more we push
        return 1.0 + (diff * (codonGC - 0.5) * 5.0);
    }

    mostFrequent(protein, targetGC) {
        let dna = '';
        let prev = null;
        let gcCount = 0;
        let nucCount = 0;

        for (const aa of protein) {
            const codons = this.codonTable.getCodonsForAA(aa);
            if (codons.length === 0) throw new Error(`Unknown amino acid: ${aa}`);
            
            const currentGC = nucCount > 0 ? gcCount / nucCount : targetGC;
            
            let best = codons[0].codon;
            let bestScore = -Infinity;

            for (const c of codons) {
                const pb = this.codonTable.getPairBias(prev, c.codon) || 1.0;
                const gw = this.getGCWeight(c.codon, currentGC, targetGC);
                
                // Score combines Frequency, Pair Bias, and GC Weight
                const score = c.freq * pb * gw;
                
                if (score > bestScore) {
                    bestScore = score;
                    best = c.codon;
                }
            }
            dna += best;
            gcCount += (best.match(/[GC]/gi) || []).length;
            nucCount += 3;
            prev = best;
        }
        return dna;
    }

    weightedRandom(protein, targetGC) {
        let dna = '';
        let prev = null;
        let gcCount = 0;
        let nucCount = 0;

        for (const aa of protein) {
            const codons = this.codonTable.getCodonsForAA(aa);
            if (codons.length === 0) throw new Error(`Unknown amino acid: ${aa}`);
            
            const currentGC = nucCount > 0 ? gcCount / nucCount : targetGC;

            const weights = codons.map(c => {
                const pb = this.codonTable.getPairBias(prev, c.codon) || 1.0;
                const gw = this.getGCWeight(c.codon, currentGC, targetGC);
                return c.freq * pb * gw;
            });

            const total = weights.reduce((s, w) => s + w, 0);
            let r = Math.random() * total;
            let chosen = codons[0].codon;
            
            for (let i = 0; i < codons.length; i++) {
                r -= weights[i];
                if (r <= 0) {
                    chosen = codons[i].codon;
                    break;
                }
            }
            dna += chosen;
            gcCount += (chosen.match(/[GC]/gi) || []).length;
            nucCount += 3;
            prev = chosen;
        }
        return dna;
    }

    balanced(protein, targetGC) {
        const usage = {};
        let dna = '';
        let prev = null;
        let gcCount = 0;
        let nucCount = 0;

        for (const aa of protein) {
            const codons = this.codonTable.getCodonsForAA(aa);
            if (codons.length === 0) throw new Error(`Unknown amino acid: ${aa}`);

            if (!usage[aa]) {
                usage[aa] = codons.map(c => ({ codon: c.codon, count: 0, targetFreq: c.freq }));
            }

            const currentGC = nucCount > 0 ? gcCount / nucCount : targetGC;
            const totalUsage = usage[aa].reduce((sum, c) => sum + c.count, 0) + 1;
            
            let bestCodon = null;
            let minDiff = Infinity;

            for (const c of usage[aa]) {
                const pb = this.codonTable.getPairBias(prev, c.codon) || 1.0;
                const gw = this.getGCWeight(c.codon, currentGC, targetGC);
                
                const currentFreq = (c.count + 1) / totalUsage;
                // Target frequency is per thousand, scaled by bias and GC weight
                const adjustedTarget = (c.targetFreq * pb * gw) / 1000;
                const diff = Math.abs(currentFreq - adjustedTarget);

                if (diff < minDiff) {
                    minDiff = diff;
                    bestCodon = c;
                }
            }

            bestCodon.count++;
            dna += bestCodon.codon;
            gcCount += (bestCodon.codon.match(/[GC]/gi) || []).length;
            nucCount += 3;
            prev = bestCodon.codon;
        }
        return dna;
    }
}
