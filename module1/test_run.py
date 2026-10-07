import sys
import main

print("Running pipeline...", flush=True)
seq = "MAKLSTDELLMAARQNLVKTPRAATVLSAPQATLVAPQATLVAPQATLVAP"
result = main.run_module1(seq)
print("\nPipeline complete!", flush=True)
