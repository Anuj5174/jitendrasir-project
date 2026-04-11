// module2/js/config.js

export const CONFIG = {
    paths: {
        codonFreq: './data/human-codon-freq.json',
        codonPairBias: './data/codon-pair-bias.json',
        foldingApi: 'http://localhost:8000/api/fold'
    },
    gc: {
        defaultTarget: 0.52,
        defaultWindowSize: 50,
        tolerance: 0.05,
        maxIterations: 5
    },
    forbidden: {
        maxIterations: 10,
        patterns: [
            { name: 'Polyadenylation signal (1)', regex: /AATAAA/gi },
            { name: 'Polyadenylation signal (2)', regex: /ATTAAA/gi },
            { name: 'Poly-A run', regex: /AAAAA/gi },
            { name: 'Poly-T run', regex: /TTTTT/gi },
            { name: 'Poly-G run', regex: /GGGG/gi },
            { name: 'Poly-C run', regex: /CCCC/gi },
            { name: 'AU-rich element', regex: /ATTTA/gi },
            { name: 'BamHI site', regex: /GGATCC/gi },
            { name: 'EcoRI site', regex: /GAATTC/gi },
            { name: 'XbaI site', regex: /TCTAGA/gi },
            { name: 'Cryptic Splice Donor (AG|GT->GU)', regex: /AGGTAAGT/gi },
            { name: 'Cryptic Splice Acceptor (CAG)', regex: /CAG/gi },
            { name: 'IRES-like motif (proxy)', regex: /CCGCC/gi },
            { name: 'GC-rich run (hairpin initiator)', regex: /G{4,}|C{4,}|(GC){4,}/gi }
        ]
    }
};
