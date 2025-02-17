"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogue
A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
This program uses several heuristics:
 - The molecule must be a valid RDKit molecule.
 - Its molecular weight should be within roughly 80–350 Da.
 - It should contain at least 2 nitrogen atoms.
 - It should be small in terms of heavy atoms (≤35 heavy atoms).
 - It should be nearly rigid (≤3 rotatable bonds) to exclude flexible nucleosides.
 - The molecule must have 1 or 2 rings (nucleobases typically are mono- or bicyclic).
 - There must be at least one aromatic heterocyclic ring (of size 5–7 members with ≥1 heteroatom).
 - In addition, the ratio of carbon atoms to nitrogen atoms must be ≤4.0
These extra filters are intended to remove false positives (e.g. molecules with excessive carbon content or extra ring systems).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    Heuristics used:
      1. Valid RDKit molecule (with sanitization).
      2. Molecular weight between 80 and 350 Da.
      3. Contains at least 2 nitrogen atoms.
      4. Contains at most 35 heavy atoms.
      5. Contains ≤3 rotatable bonds (reflecting rigidity).
      6. Contains 1 or 2 rings (typical for nucleobases).
      7. Contains at least one aromatic heterocyclic ring of size 5–7 (with ≥1 heteroatom).
      8. The ratio of carbon atoms to nitrogen atoms is ≤4.0.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a nucleobase analogue, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Check molecular weight (80-350 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside expected range (80–350 Da)"
    
    # Count nitrogen atoms; must have at least 2.
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogen < 2:
        return False, "Too few nitrogen atoms; a nucleobase analogue should contain at least 2 nitrogen atoms"
    
    # Count heavy atoms; require at most 35.
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 35:
        return False, f"Too many heavy atoms ({heavy_atom_count}); nucleobase analogues are typically small (≤35 heavy atoms)"
    
    # Check rotatable bonds; should be ≤3, as nucleobases are rigid.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 3:
        return False, f"Too many rotatable bonds ({n_rotatable}); nucleobase analogues are expected to be rigid (≤3 rotatable bonds)"
    
    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    n_rings = len(atom_rings)
    if n_rings < 1:
        return False, "No ring structures found; nucleobase analogues are heterocyclic"
    if n_rings > 2:
        return False, f"Too many rings ({n_rings}); natural nucleobases are usually mono- or bicyclic"
    
    # Look for at least one appropriate aromatic heterocyclic ring.
    heterocycle_found = False
    for ring in atom_rings:
        ring_size = len(ring)
        if ring_size in (5, 6, 7):
            # Check that every atom in the ring is aromatic.
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                # Count heteroatoms (non-carbon atoms) in the ring.
                hetero_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6)
                if hetero_count >= 1:
                    heterocycle_found = True
                    break
    if not heterocycle_found:
        return False, "No appropriate aromatic heterocyclic ring (size 5–7 with at least one heteroatom) found"
    
    # Count carbon atoms.
    n_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Check carbon-to-nitrogen ratio.
    if n_nitrogen > 0 and (n_carbon / n_nitrogen) > 4.0:
        return False, f"High carbon-to-nitrogen ratio ({n_carbon}/{n_nitrogen} = {n_carbon/n_nitrogen:.2f}); nucleobase analogues are generally nitrogen rich"
    
    return True, ("Molecule exhibits features of a nucleobase analogue: appropriate molecular weight, nitrogen count, small size, rigidity, "
                  "limited ring system, an aromatic heterocyclic ring (5–7 members), and a balanced C/N ratio.")

# Example test calls (uncomment below to test)
# if __name__ == "__main__":
#     test_smiles = [
#         "C=12C(=NC=NC1NCC(=O)N)NC=N2",  # N(6)-carbamoylmethyladenine (True)
#         "Nc1nc(=O)[nH]cc1CO",            # 5-(hydroxymethyl)cytosine (True)
#         "COC(=O)Nc1nc2cc(ccc2[nH]1)C(=O)c1ccccc1"  # mebendazole (False)
#     ]
#     for sm in test_smiles:
#         res, reason = is_nucleobase_analogue(sm)
#         print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")