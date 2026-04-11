// module2/js/frame-builder.js

export class FrameBuilder {
    static KOZAK = 'GCCACC';

    build(cds, stopCodon = 'TGA') {
        const startCodon = 'ATG';
        // Input CDS might already have ATG, let's normalize
        let normalizedCds = cds.toUpperCase();
        if (normalizedCds.startsWith('ATG')) {
            normalizedCds = normalizedCds.substring(3);
        }
        
        // Remove trailing stop if present to ensure we use the requested one
        if (normalizedCds.endsWith('TGA') || normalizedCds.endsWith('TAA') || normalizedCds.endsWith('TAG')) {
            normalizedCds = normalizedCds.substring(0, normalizedCds.length - 3);
        }

        const fullSequence = FrameBuilder.KOZAK + startCodon + normalizedCds + stopCodon;
        
        return {
            kozak: FrameBuilder.KOZAK,
            startCodon,
            cds: normalizedCds,
            stopCodon,
            fullSequence
        };
    }
}
