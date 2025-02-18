"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: 16β–hydroxy steroid
Definition: A 16–hydroxy steroid in which the hydroxy group at position 16 has a beta–configuration.
This program uses a simple heuristic: first it checks for a steroid–like fused four–ring system 
(with rings of size 5 or 6) and then looks for a chiral carbon in a five–membered ring (usually ring D)
bearing an –OH. Finally, it heuristically assumes that if the chiral tag for that carbon is CHI_TETRAHEDRAL_CW, 
then the hydroxy is beta–oriented.
Note: This method is only a heuristic and may fail for many edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16β–hydroxy steroid based on its SMILES string.
    
    A valid 16β–hydroxy steroid is defined (heuristically) as a molecule with a steroid nucleus 
    (a fused ring system of 4 rings predominantly of size 5 or 6) and a chiral carbon in one of the 
    five-membered rings (typically the D–ring) which is bonded to a hydroxyl (-OH) group. 
    Furthermore, we assume that if the chiral tag of that carbon is CHI_TETRAHEDRAL_CW the –OH is in beta–configuration.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a 16β–hydroxy steroid, False otherwise.
        str: Reason explaining the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Heuristic: Assume a steroid nucleus should contain at least 4 rings of size 5 or 6.
    steroid_rings = [r for r in rings if len(r) in (5, 6)]
    if len(steroid_rings) < 4:
        return False, f"Found only {len(steroid_rings)} rings of size 5 or 6; steroid nucleus typically has 4 fused rings."
    
    # Focus on five-membered rings (commonly the D-ring in steroids)
    five_membered_rings = [r for r in rings if len(r) == 5]
    if not five_membered_rings:
        return False, "No five-membered ring found; expected D-ring in steroid nucleus."
    
    # Look for a chiral carbon in a five-membered ring that bears an –OH.
    candidate_found = False
    hydroxyl_atom = None
    for ring in five_membered_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Consider only carbon atoms
            if atom.GetAtomicNum() != 6:
                continue
            # Check if the atom is marked as chiral.
            if atom.GetChiralTag() == rdchem.ChiralType.CHI_UNSPECIFIED:
                continue
            # Look among its neighbors to see if any is a hydroxyl oxygen (with a single bond)
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.SINGLE:
                        candidate_found = True
                        hydroxyl_atom = atom
                        break
            if candidate_found:
                break
        if candidate_found:
            break

    if not candidate_found:
        return False, "No chiral carbon in a five-membered ring with an attached hydroxyl group was detected."
    
    # Attempt to judge the configuration.
    # NOTE: RDKit assigns chiral tags but does not directly label beta vs. alpha.
    # Here we use a heuristic: assume that if the chiral tag of the hydroxyl-bearing carbon 
    # is CHI_TETRAHEDRAL_CW then it is beta; otherwise, we assume it is not.
    if hydroxyl_atom.GetChiralTag() != rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
        return False, "The chiral center bearing the hydroxyl group does not have the expected (beta) configuration (heuristic)."
    
    return True, "Molecule contains a steroid nucleus with a 16-hydroxyl group in beta configuration (heuristic)."

# Example usage (uncomment to test):
# smiles_example = "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"  # 16β-hydroxytestosterone
# print(is_16beta_hydroxy_steroid(smiles_example))