"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: Wax ester, defined as a fatty acid ester resulting from 
the condensation of the carboxy group of a fatty acid with the alcoholic hydroxy group 
of a fatty alcohol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is formed by the condensation of a fatty acid with a fatty alcohol,
    resulting in a single ester bond (C(=O)O) linking two long aliphatic chains.
    
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

    # Define an ester group SMARTS pattern.
    # The pattern "[CX3](=O)[OX2]" matches three atoms:
    #  1. A carbon (the carbonyl carbon)
    #  2. Its double-bonded oxygen atom
    #  3. The oxygen atom that links to the fatty alcohol fragment.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester matches; expected exactly 1 for a wax ester."
        
    # Unpack the three indices from the match.
    # ester_matches[0] contains (carbonyl_idx, dbl_bonded_oxygen_idx, alcohol_oxygen_idx)
    carbonyl_idx, dbl_ox_idx, alcohol_ox_idx = ester_matches[0]
    
    # Get the bond between the carbonyl carbon and the alcohol oxygen.
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, alcohol_ox_idx)
    if bond is None:
        return False, "No bond found between carbonyl carbon and ester oxygen."
    bond_idx = bond.GetIdx()
    
    # Fragment the molecule on the ester bond.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    
    if len(frags) != 2:
        return False, f"Fragmentation yielded {len(frags)} fragments; expected 2 fragments."
    
    # Helper function: count the number of carbon atoms in a fragment ignoring dummy atoms.
    def count_carbons(fragment):
        return sum(1 for atom in fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Identify the fatty acid (acyl) and fatty alcohol fragments.
    # A fatty acid fragment should contain a carbonyl group, so we search for "[CX3](=O)".
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)")
    frag_acyl = None
    frag_alcohol = None

    for frag in frags:
        # Check if the fragment has a carbonyl group.
        if frag.HasSubstructMatch(acid_pattern):
            frag_acyl = frag
        else:
            frag_alcohol = frag
            
    # If our assignment did not work via pattern matching, assign based on dummy atoms.
    if frag_acyl is None or frag_alcohol is None:
        for frag in frags:
            # Look for dummy atoms (atomic number 0) which should be on the alcohol side.
            if any(atom.GetAtomicNum() == 0 for atom in frag.GetAtoms()):
                frag_alcohol = frag
            else:
                frag_acyl = frag

    if frag_acyl is None or frag_alcohol is None:
        return False, "Could not identify both fatty acid and fatty alcohol fragments."
        
    # Count carbons in each fragment.
    acid_carbons = count_carbons(frag_acyl)
    alcohol_carbons = count_carbons(frag_alcohol)
    
    # Set a threshold for a long aliphatic chain (typically at least 6 carbons are required).
    min_chain_length = 6
    if acid_carbons < min_chain_length:
        return False, f"Fatty acid fragment too short ({acid_carbons} carbons; minimum {min_chain_length} required)."
    if alcohol_carbons < min_chain_length:
        return False, f"Fatty alcohol fragment too short ({alcohol_carbons} carbons; minimum {min_chain_length} required)."
        
    # Optionally, check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a typical wax ester."
    
    return True, f"Single ester group found with fatty acid fragment ({acid_carbons} C) and fatty alcohol fragment ({alcohol_carbons} C)."
    
# Example usage:
# result, reason = is_wax_ester("CCCCCCCCCCCCCCCC(=O)OCCCCCCCCC")
# print(result, reason)