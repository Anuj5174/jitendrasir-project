import { CONFIG } from './config.js';

export class CodonTable {
    constructor() {
        this.data = null;
        this.aaToCodons = {};
        this.codonToAA = {};
    }

    async init() {
        const response = await fetch(CONFIG.paths.codonFreq);
        this.data = await response.json();

        // Try to load codon-pair bias matrix (optional)
        try {
            const cpbResp = await fetch(CONFIG.paths.codonPairBias);
            this.pairBias = await cpbResp.json();
        } catch (e) {
            this.pairBias = { default: 1.0 };
        }

        for (const [codon, info] of Object.entries(this.data)) {
            const aa = info.aa;
            if (!this.aaToCodons[aa]) {
                this.aaToCodons[aa] = [];
            }
            this.aaToCodons[aa].push({
                codon,
                freq: info.freq
            });
            this.codonToAA[codon] = aa;
        }

        // Sort codons by frequency for each AA
        for (const aa in this.aaToCodons) {
            this.aaToCodons[aa].sort((a, b) => b.freq - a.freq);
        }
    }

    getCodonsForAA(aa) {
        return this.aaToCodons[aa] || [];
    }

    getMostFrequent(aa) {
        const codons = this.getCodonsForAA(aa);
        return codons.length > 0 ? codons[0].codon : null;
    }

    getRelativeAdaptiveness(codon) {
        const aa = this.codonToAA[codon];
        if (!aa) return 0;
        const maxFreq = this.aaToCodons[aa][0].freq;
        return this.data[codon].freq / maxFreq;
    }

    getAAForCodon(codon) {
        return this.codonToAA[codon];
    }

    getPairBias(prevCodon, nextCodon) {
        if (!this.pairBias) return 1.0;
        if (!prevCodon || !nextCodon) return this.pairBias.default || 1.0;
        const key = `${prevCodon}-${nextCodon}`;
        return this.pairBias[key] || this.pairBias.default || 1.0;
    }
}
