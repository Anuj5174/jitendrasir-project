// module2/js/cai-calculator.js

export class CAICalculator {
    constructor(codonTable) {
        this.codonTable = codonTable;
    }

    calculate(dna) {
        if (!dna || dna.length === 0) return 0;
        
        const codons = [];
        for (let i = 0; i < dna.length; i += 3) {
            codons.push(dna.substring(i, i + 3));
        }

        let logSum = 0;
        let count = 0;

        for (const codon of codons) {
            const w = this.codonTable.getRelativeAdaptiveness(codon);
            if (w > 0) {
                logSum += Math.log(w);
                count++;
            }
        }

        if (count === 0) return 0;
        
        // CAI is the geometric mean of w values
        // geometric mean = exp( (1/L) * sum(ln(w_i)) )
        const cai = Math.exp(logSum / count);
        return parseFloat(cai.toFixed(3));
    }
}
