// module2/js/verify.js

export class Verifier {
    constructor(codonTable) {
        this.codonTable = codonTable;
    }

    verify(dna, originalProtein) {
        // Remove Kozak and stop codon for verification
        // Assumes GCCACC ATG ... STOP
        let cds = dna.toUpperCase();
        if (cds.startsWith('GCCACC')) {
            cds = cds.substring(6);
        }
        
        let expectedProtein = originalProtein;
        if (!expectedProtein.startsWith('M')) {
            expectedProtein = 'M' + expectedProtein;
        }

        let translated = '';
        const mismatches = [];

        for (let i = 0; i < cds.length; i += 3) {
            const codon = cds.substring(i, i + 3);
            if (codon.length < 3) break; // Trailing bases

            const aa = this.codonTable.getAAForCodon(codon);
            if (aa === '*') {
                // Stop codon reached
                break;
            }
            
            translated += aa;
            
            const proteinIdx = i / 3;
            if (aa !== expectedProtein[proteinIdx]) {
                mismatches.push({
                    pos: proteinIdx + 1,
                    expected: expectedProtein[proteinIdx],
                    actual: aa,
                    codon: codon
                });
            }
        }

        const success = translated === expectedProtein;
        
        return {
            success,
            translated,
            mismatches,
            lengthMatch: translated.length === expectedProtein.length
        };
    }
}
