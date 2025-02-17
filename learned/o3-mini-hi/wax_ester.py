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
    which results in a single ester bond (C(=O)O) linking two long aliphatic chains.
    
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
    # This pattern looks for a carbonyl carbon (with a double bond O) attached to an oxygen.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups; expected exactly 1 for a wax ester."
        
    # Get the match indices.
    # In the SMARTS "[CX3](=O)[OX2]", the match returns two atoms:
    #   match[0] --> the carbonyl carbon
    #   match[1] --> the ester oxygen (that links to the fatty alcohol fragment)
    carbonyl_idx, oxy_idx = ester_matches[0]
    
    # Identify the two fragments: one containing the carbonyl carbon (the fatty acid part)
    # and the other containing the ester oxygen (the fatty alcohol part).
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, oxy_idx)
    if bond is None:
        return False, "No bond found between carbonyl carbon and ester oxygen."
    bond_idx = bond.GetIdx()
    
    # Break the ester bond by fragmenting the molecule on that bond.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    
    if len(frags) != 2:
        return False, f"Fragmentation yielded {len(frags)} fragments; expected 2 fragments."
    
    # Helper function: count the number of carbon atoms in a fragment ignoring dummy atoms.
    def count_carbons(fragment):
        return sum(1 for atom in fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Determine which fragment is which:
    # We check for the presence of our key atoms. In the fragmented molecule, dummy atoms (atomic number 0)
    # appear at the broken bond. We use the original indices to figure out which fragment contains the carbonyl carbon.
    frag_acyl = None
    frag_alcohol = None
    for frag in frags:
        atom_ids = [atom.GetIdx() for atom in frag.GetAtoms() if atom.GetAtomicNum()!=0]
        # Get SMILES for debugging; alternatively, check if the fragment has a carbonyl functionality.
        frag_smiles = Chem.MolToSmiles(frag)
        # If the fragment contains a carbonyl pattern, assign it as the fatty acid part.
        # We look for a "[CX3](=O)" pattern.
        acid_pattern = Chem.MolFromSmarts("[CX3](=O)")
        if frag.HasSubstructMatch(acid_pattern):
            frag_acyl = frag
        else:
            frag_alcohol = frag
            
    # If our assignment did not work via patterns, assign based on which fragment has a dummy atom 
    # attached to an oxygen (typical of the alcohol side) and which one has a carbonyl.
    if frag_acyl is None or frag_alcohol is None:
        # Fall back: assign based on the ester oxygen presence.
        for frag in frags:
            frag_atoms = frag.GetAtoms()
            has_ester_oxygen = any(atom.GetAtomicNum() == 8 and atom.GetSymbol() == "O" and atom.GetDegree() > 0 
                                     for atom in frag_atoms)
            if has_ester_oxygen and frag_alcohol is None:
                frag_alcohol = frag
            else:
                frag_acyl = frag
                
    if frag_acyl is None or frag_alcohol is None:
        return False, "Could not identify both fatty acid and fatty alcohol fragments."
        
    # Count carbons in each fragment.
    acid_carbons = count_carbons(frag_acyl)
    alcohol_carbons = count_carbons(frag_alcohol)
    
    # Set a threshold for a long aliphatic chain (typically at least 6 carbons).
    min_chain_length = 6
    if acid_carbons < min_chain_length:
        return False, f"Fatty acid fragment too short ({acid_carbons} carbons; minimum {min_chain_length} required)."
    if alcohol_carbons < min_chain_length:
        return False, f"Fatty alcohol fragment too short ({alcohol_carbons} carbons; minimum {min_chain_length} required)."
        
    # Optionally, one could add further checks on overall molecular weight or lipophilicity.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a typical wax ester."
    
    return True, f"Single ester group found with fatty acid fragment ({acid_carbons} C) and fatty alcohol fragment ({alcohol_carbons} C)."
    
# Example usage:
# result, reason = is_wax_ester("CCCCCCCCCCCCCCCC(=O)OCCCCCCCCC")
# print(result, reason)