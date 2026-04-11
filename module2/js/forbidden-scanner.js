import { CONFIG } from './config.js';

export class ForbiddenScanner {
    constructor(codonTable) {
        this.codonTable = codonTable;
        this.patterns = CONFIG.forbidden.patterns;
    }

    scan(dna) {
        const matches = [];
        for (const p of this.patterns) {
            let match;
            while ((match = p.regex.exec(dna)) !== null) {
                matches.push({
                    name: p.name,
                    index: match.index,
                    sequence: match[0]
                });
            }
            p.regex.lastIndex = 0; // Reset regex index
        }
        return matches;
    }

    fix(dna) {
        let currentDna = dna;
        const codons = [];
        for (let i = 0; i < currentDna.length; i += 3) {
            codons.push(currentDna.substring(i, i + 3));
        }

        let iterations = 0;
        const maxIterations = 10;

        while (iterations < maxIterations) {
            const matches = this.scan(currentDna);
            if (matches.length === 0) break;

            for (const match of matches) {
                // Find all codons that overlap with this match
                const startCodonIdx = Math.floor(match.index / 3);
                const endCodonIdx = Math.ceil((match.index + match.sequence.length) / 3);

                for (let i = startCodonIdx; i < endCodonIdx; i++) {
                    const aa = this.codonTable.getAAForCodon(codons[i]);
                    const synonymous = this.codonTable.getCodonsForAA(aa);

                    if (synonymous.length > 1) {
                        for (const s of synonymous) {
                            if (s.codon === codons[i]) continue;

                            const testCodons = [...codons];
                            testCodons[i] = s.codon;
                            const testDna = testCodons.join('');
                            
                            // Check if this fix resolved at least one occurrence of THIS pattern
                            const newMatches = this.scan(testDna).filter(m => m.name === match.name);
                            if (newMatches.length < matches.filter(m => m.name === match.name).length) {
                                codons[i] = s.codon;
                                currentDna = testDna;
                                break;
                            }
                        }
                    }
                }
            }
            iterations++;
        }

        return currentDna;
    }
}
