// module2/js/cpg-audit.js

export class CpGAudit {
    audit(dna) {
        if (!dna) return { count: 0, density: 0, flag: false };
        
        const matches = dna.match(/CG/gi) || [];
        const count = matches.length;
        const density = (count * 2 / dna.length) * 100;
        const flag = density > 3.0;

        return {
            count,
            density: parseFloat(density.toFixed(2)),
            flag
        };
    }
}
