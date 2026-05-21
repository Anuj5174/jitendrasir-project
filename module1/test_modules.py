import sys
sys.path.insert(0, 'd:/ONEHEALTHONENATION/jitendrasir-project/module1')

print("--- Import checks ---")
import antigenicity; print("antigenicity    OK")
import protparam;    print("protparam       OK")
import homology;     print("homology        OK")
import immune_export; print("immune_export   OK")
import fusion;       print("fusion          OK")

print("\n--- VaxiJen-like antigenicity test ---")
from antigenicity import predict_antigenicity
from default_config import DEFAULT_CONFIG
r = predict_antigenicity('KVGAAFTLI', DEFAULT_CONFIG)
print("Score:", r["score"], "  is_antigen:", r["is_antigen"])
print("Note:", r["note"])

print("\n--- ProtParam test ---")
from protparam import analyze_construct
test_seq = 'MAKLSTDELLDAFKEMTLLELSDFVKK'
p = analyze_construct(test_seq)
print("MW:", p["molecular_weight_kda"], "kDa")
print("pI:", p["theoretical_pi"])
print("II:", p["instability_index"]["value"], "->", p["instability_index"]["classification"])
print("GRAVY:", p["gravy"])

print("\n--- Fusion test (full construct) ---")
from fusion import fuse_epitopes
epitopes = [
    {"peptide": "KVGAAFTLI",  "type": "MHC-I"},
    {"peptide": "YPLLWSFAMG", "type": "MHC-I"},
    {"peptide": "LLNLRSRLA",  "type": "MHC-II"},
    {"peptide": "FAMGVATTI",  "type": "B-cell"},
]
construct = fuse_epitopes(epitopes, DEFAULT_CONFIG)
print("Length:", len(construct), "aa")
print("Starts with:", construct[:20])
print("Ends with:  ", construct[-10:])

print("\n--- Immune report test ---")
from immune_export import generate_immune_report
report = generate_immune_report(construct, epitopes, DEFAULT_CONFIG)
ig = report["heuristic_immune_response"]["interpretation"]
print("IgM:", ig["IgM"], " IgG:", ig["IgG"], " CD8:", ig["CD8"], " Memory:", ig["Memory"])
print("Readiness:", repr(report["readiness_status"]))

print("\nAll tests passed!")
