"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax ester, defined as a fatty acid ester resulting from 
the condensation of the carboxy group of a fatty acid with the alcoholic hydroxy group 
of a fatty alcohol.
Improved version: rejects molecules with additional ester groups and ensures that 
both fragments are long unfunctionalized aliphatic chains.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is produced by the condensation of a fatty acid (which supplies the carbonyl group)
    with a fatty alcohol, forming a single ester bond. Both resulting fragments should have long 
    unfunctionalized (aliphatic) chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a wax ester, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a single ester bond using SMARTS.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Reject if not exactly one ester group is found.
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester matches; expected exactly 1 for a wax ester."
    
    # Unpack the match: (carbonyl carbon, double-bonded oxygen, and linking oxygen)
    carbonyl_idx, dbl_ox_idx, alcohol_ox_idx = ester_matches[0]
    
    # Verify the bond between the carbonyl carbon and the alcohol oxygen exists.
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, alcohol_ox_idx)
    if bond is None:
        return False, "No bond found between carbonyl carbon and ester oxygen."
    bond_idx = bond.GetIdx()
    
    # Fragment the molecule by breaking the ester bond.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    if len(frags) != 2:
        return False, f"Fragmentation yielded {len(frags)} fragments; expected 2 fragments."
    
    # Define a helper to count aliphatic (non-dummy, non-aromatic) carbons.
    def count_aliphatic_carbons(fragment):
        count = 0
        for atom in fragment.GetAtoms():
            # Exclude dummy atoms (atomic number 0).
            if atom.GetAtomicNum() != 6:
                continue
            # Exclude aromatic carbons.
            if atom.GetIsAromatic():
                continue
            count += 1
        return count
    
    # Define a SMARTS to detect a carbonyl (for identifying the acid fragment).
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)")
    
    frag_acyl = None
    frag_alcohol = None
    for frag in frags:
        if frag.HasSubstructMatch(acid_pattern):
            frag_acyl = frag
        else:
            frag_alcohol = frag
    # If assignment failed, try assigning by dummy atoms:
    if frag_acyl is None or frag_alcohol is None:
        for frag in frags:
            # The fatty alcohol fragment should be the one that has a dummy atom (due to fragmentation)
            if any(atom.GetAtomicNum() == 0 for atom in frag.GetAtoms()):
                frag_alcohol = frag
            else:
                frag_acyl = frag
                
    if frag_acyl is None or frag_alcohol is None:
        return False, "Could not identify both fatty acid and fatty alcohol fragments."
    
    acid_carbons = count_aliphatic_carbons(frag_acyl)
    alcohol_carbons = count_aliphatic_carbons(frag_alcohol)
    
    # Define minimum thresholds for the two fragments 
    # (typical wax esters have a long fatty acid chain and a long fatty alcohol chain).
    min_acid_carbons = 12
    min_alcohol_carbons = 8
    
    if acid_carbons < min_acid_carbons:
        return False, f"Fatty acid fragment too short ({acid_carbons} aliphatic carbons; minimum {min_acid_carbons} required)."
    if alcohol_carbons < min_alcohol_carbons:
        return False, f"Fatty alcohol fragment too short ({alcohol_carbons} aliphatic carbons; minimum {min_alcohol_carbons} required)."
    
    # Optionally, check that the fragments do not contain aromatic atoms (which would indicate non-aliphatic groups).
    def contains_aromatic(fragment):
        return any(atom.GetAtomicNum() == 6 and atom.GetIsAromatic() for atom in fragment.GetAtoms())
    
    if contains_aromatic(frag_acyl):
        return False, "Fatty acid fragment contains aromatic groups."
    if contains_aromatic(frag_alcohol):
        return False, "Fatty alcohol fragment contains aromatic groups."
    
    # Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a typical wax ester."
    
    return True, f"Single ester group found with fatty acid fragment ({acid_carbons} aliphatic C) and fatty alcohol fragment ({alcohol_carbons} aliphatic C)."
    
# Example usage:
# result, reason = is_wax_ester("CCCCCCCCCCCCCCCC(=O)OCCCCCCCCC")
# print(result, reason)