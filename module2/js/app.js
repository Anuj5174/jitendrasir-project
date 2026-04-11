// module2/js/app.js

import { CONFIG } from './config.js';
import { CodonTable } from './codon-table.js';
import { ReverseTranslator } from './reverse-translate.js';
import { GCBalancer } from './gc-balancer.js';
import { ForbiddenScanner } from './forbidden-scanner.js';
import { CpGAudit } from './cpg-audit.js';
import { CAICalculator } from './cai-calculator.js';
import { FrameBuilder } from './frame-builder.js';
import { Verifier } from './verify.js';
import { Exporter } from './exporter.js';

class App {
    constructor() {
        this.codonTable = new CodonTable();
        this.reverseTranslator = null;
        this.gcBalancer = null;
        this.forbiddenScanner = null;
        this.cpgAudit = new CpGAudit();
        this.caiCalculator = null;
        this.frameBuilder = new FrameBuilder();
        this.verifier = null;
        this.exporter = new Exporter();

        this.initialized = false;
        this.currentData = null;

        this.bindEvents();
    }

    async init() {
        if (this.initialized) return;
        await this.codonTable.init();
        
        this.reverseTranslator = new ReverseTranslator(this.codonTable);
        this.gcBalancer = new GCBalancer(this.codonTable);
        this.forbiddenScanner = new ForbiddenScanner(this.codonTable);
        this.caiCalculator = new CAICalculator(this.codonTable);
        this.verifier = new Verifier(this.codonTable);
        
        this.initialized = true;
    }

    bindEvents() {
        document.getElementById('generateBtn').addEventListener('click', () => this.runPipeline());
        document.getElementById('copyBtn').addEventListener('click', () => this.copySequence());
        document.getElementById('exportBtn').addEventListener('click', () => this.exportJson());
    }

    async runPipeline() {
        await this.init();
        
        let proteinRaw = document.getElementById('proteinInput').value.trim().toUpperCase();
        const protein = proteinRaw.replace(/[^A-Z]/g, ''); // Removes spaces and special characters
        if (!protein) {
            alert('Please enter a protein sequence.');
            return;
        }

        const strategy = document.getElementById('strategySelect').value;
        const targetGC = parseInt(document.getElementById('targetGC').value) / 100;
        const stopCodon = document.getElementById('stopCodon').value;

        document.getElementById('progressSection').classList.remove('hidden');
        this.resetUI();

        try {
            let dna = '';
            let cai = 0;
            let gc = 0;
            let bestConstruct = null;
            let bestMfe = Infinity;
            let finalCai = 0;
            let finalGc = 0;

            const CANDIDATES_COUNT = CONFIG.optimization.candidatesCount;
            console.log(`Starting Research-Grade optimization (Candidates: ${CANDIDATES_COUNT})...`);

            const candidates = [];
            for (let c = 0; c < CANDIDATES_COUNT; c++) {
                console.log(`Generating candidate ${c + 1}...`);
                let dna = this.reverseTranslator.translate(protein, strategy, targetGC);
                dna = this.gcBalancer.balance(dna, targetGC);
                dna = this.forbiddenScanner.fix(dna);
                
                const construct = this.frameBuilder.build(dna, stopCodon);
                const cai = this.caiCalculator.calculate(construct.cds);
                
                // Only consider candidates meeting the Research-grade floor
                if (cai >= CONFIG.optimization.caiFloor) {
                    candidates.push({ construct, cai });
                }
            }

            if (candidates.length === 0) {
                throw new Error(`Failed to generate candidates meeting CAI ${CONFIG.optimization.caiFloor} threshold. Try a different strategy.`);
            }

            // 4. Structural Selection (Fold all and pick best)
            this.updateStep(4, 'active');
            for (const cand of candidates) {
                let mfe = 0;
                try {
                    const rnaSeq = cand.construct.cds.replace(/T/g, 'U');
                    const resp = await fetch(CONFIG.paths.foldingApi, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ sequence: rnaSeq })
                    });
                    if (resp.ok) {
                        const fold = await resp.json();
                        mfe = fold.mfe;
                        cand.fold = fold;
                    }
                } catch (e) {
                    console.warn("Structural selection fallback to heuristic...");
                }

                if (mfe < bestMfe) {
                    bestMfe = mfe;
                    bestConstruct = cand.construct;
                    finalCai = cand.cai;
                    bestFold = cand.fold;
                }
            }
            
            finalGc = this.gcBalancer.calculateGC(bestConstruct.cds);
            const cpg = this.cpgAudit.audit(bestConstruct.cds);
            this.updateStep(4, 'completed');

            // Folding: call local Python/ViennaRNA folding service (if available)
            let foldResult = null;
            try {
                const rnaSeq = construct.cds.replace(/T/g, 'U');
                const resp = await fetch(CONFIG.paths.foldingApi, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ sequence: rnaSeq })
                });
                if (resp.ok) {
                    foldResult = await resp.json();
                } else {
                    console.warn('Folding service returned', resp.status);
                }
            } catch (e) {
                console.warn('Folding service unreachable', e.message);
            }

            // 5. Verification
            this.updateStep(5, 'active');
            const verification = this.verifier.verify(bestConstruct.fullSequence, protein);
            const researchValidation = this.verifier.validateResearchTargets({
                cai: finalCai,
                gc: finalGc,
                mfe: bestMfe
            }, bestConstruct.cds.length);
            this.updateStep(5, 'completed');

            // Update Metrics
            this.updateMetrics(finalCai, finalGc, cpg, verification, bestFold);
            this.displaySequence(bestConstruct.fullSequence);
            
            this.currentData = {
                protein,
                strategy,
                metrics: { cai: finalCai, gc_percent: (finalGc * 100).toFixed(1), cpg },
                construct: bestConstruct,
                verification,
                researchValidation
            };

        } catch (error) {
            console.error(error);
            alert(`Pipeline Error: ${error.message}`);
        }
    }

    resetUI() {
        ['step1', 'step2', 'step3', 'step4', 'step5'].forEach(id => {
            document.getElementById(id).className = 'step';
        });
        document.getElementById('sequenceDisplay').innerHTML = '';
        document.getElementById('verifyStatus').className = 'status-badge';
        document.getElementById('verifyStatus').innerText = 'PENDING';
    }

    updateStep(num, status) {
        document.getElementById(`step${num}`).className = `step ${status}`;
    }

    updateMetrics(cai, gc, cpg, verification) {
        document.getElementById('caiValue').innerText = cai.toFixed(3);
        document.getElementById('gcValue').innerText = `${(gc * 100).toFixed(1)}%`;
        document.getElementById('cpgValue').innerText = `${cpg.density}%`;
        // optional folding result
        const mfeEl = document.getElementById('mfeValue');
        mfeEl.innerText = (arguments[4] && arguments[4].mfe !== undefined) ? `${arguments[4].mfe} kcal/mol` : 'N/A';
        
        const vBadge = document.getElementById('verifyStatus');
        if (verification.success) {
            vBadge.innerText = 'PASSED';
            vBadge.className = 'status-badge status-pass';
        } else {
            vBadge.innerText = 'FAILED';
            vBadge.className = 'status-badge status-fail';
        }
    }

    displaySequence(seq) {
        const container = document.getElementById('sequenceDisplay');
        container.innerHTML = `>optimized_cds\n${seq}`;
        container.style.fontFamily = "'Fira Code', monospace";
        container.style.whiteSpace = 'pre-wrap';
        container.style.wordBreak = 'break-all';
        container.style.color = '#e2e8f0';
        container.style.padding = '1rem';
        
        console.log("--- MODULE 2 OUTPUT ---");
        console.log(`>optimized_cds\n${seq}`);
    }

    copySequence() {
        if (!this.currentData) return;
        navigator.clipboard.writeText(this.currentData.construct.fullSequence);
        alert('Sequence copied to clipboard!');
    }

    exportJson() {
        if (!this.currentData) return;
        this.exporter.export(this.currentData);
    }
}

new App();
