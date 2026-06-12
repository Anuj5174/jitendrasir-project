import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg') # Headless
import matplotlib.pyplot as plt

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import PPBuilder
except ImportError:
    pass

def _is_in_favored(phi, psi):
    """Very simplified geometric heuristic for Favored vs Outlier for Ramachandran plot"""
    if phi is None or psi is None: return False
    # Convert to degrees for easier logic
    phi_d = math.degrees(phi)
    psi_d = math.degrees(psi)
    
    # Beta sheet region (top left)
    if -150 <= phi_d <= -50 and 90 <= psi_d <= 180: return True
    # Alpha helix (bottom left)
    if -140 <= phi_d <= -40 and -70 <= psi_d <= -30: return True
    # Left handed helix (top right)
    if 40 <= phi_d <= 80 and 20 <= psi_d <= 80: return True
    
    return False

def _is_in_allowed(phi, psi):
    """Broader allowed bounds"""
    if phi is None or psi is None: return False
    phi_d = math.degrees(phi)
    psi_d = math.degrees(psi)
    
    # Beta
    if -180 <= phi_d <= -30 and 50 <= psi_d <= 180: return True
    if -180 <= phi_d <= -30 and -180 <= psi_d <= -150: return True # Wraparound
    # Alpha
    if -180 <= phi_d <= -30 and -100 <= psi_d <= 0: return True
    # L-helix
    if 30 <= phi_d <= 100 and -20 <= psi_d <= 100: return True
    
    return False

def generate_ramachandran_plot(pdb_file: str, output_image: str) -> dict:
    """
    Parses a PDB file to calculate phi/psi angles, generates a Ramachandran plot, 
    and returns geometry statistics (Favored, Allowed, Outlier).
    """
    if not os.path.exists(pdb_file):
        return {"status": "error", "message": "PDB file not found"}

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('vaccine', pdb_file)
    except NameError:
        return {"status": "error", "message": "Biopython not installed"}

    # Extract phi/psi angles
    ppb = PPBuilder()
    phi_psi = []
    
    for pp in ppb.build_peptides(structure):
        angles = pp.get_phi_psi_list()
        for phi, psi in angles:
            if phi is not None and psi is not None:
                phi_psi.append((phi, psi))

    if not phi_psi:
        return {"status": "error", "message": "No backbone angles found in PDB"}

    # Calculate statistics
    total = len(phi_psi)
    favored = sum(1 for p, s in phi_psi if _is_in_favored(p, s))
    allowed = sum(1 for p, s in phi_psi if not _is_in_favored(p, s) and _is_in_allowed(p, s))
    outliers = total - favored - allowed

    favored_pct = round((favored / total) * 100, 1)
    allowed_pct = round((allowed / total) * 100, 1)
    outlier_pct = round((outliers / total) * 100, 1)

    # Plotting
    plt.figure(figsize=(6, 6))
    
    # Draw simple background regions
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    
    # Extract arrays in degrees
    x = [math.degrees(phi) for phi, psi in phi_psi]
    y = [math.degrees(psi) for phi, psi in phi_psi]

    # Scatter plot
    plt.scatter(x, y, s=15, c='#3b82f6', edgecolors='#1e3a8a', alpha=0.7)
    
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.title('Ramachandran Plot (PROCHECK-style)')
    plt.xlabel('Phi ($\phi$) Angle')
    plt.ylabel('Psi ($\psi$) Angle')
    
    plt.grid(color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
    
    # Text overlay
    plt.text(-170, 160, f"Favored: {favored_pct}%", fontsize=10, color='green' if favored_pct > 90 else 'orange')
    plt.text(-170, 140, f"Allowed: {allowed_pct}%", fontsize=10)
    plt.text(-170, 120, f"Outliers: {outlier_pct}%", fontsize=10, color='red' if outlier_pct > 5 else 'black')

    plt.tight_layout()
    plt.savefig(output_image, dpi=150)
    plt.close()

    return {
        "status": "success",
        "favored": favored_pct,
        "allowed": allowed_pct,
        "outliers": outlier_pct,
        "plot_path": output_image,
        "interpretation": "Excellent model" if favored_pct > 90 else "Acceptable model" if favored_pct > 80 else "Poor model quality"
    }

if __name__ == "__main__":
    res = generate_ramachandran_plot("output/structure.pdb", "output/ramachandran.png")
    print(res)
