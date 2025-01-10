"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester
Definition: A fatty acid ester resulting from the formal condensation of the carboxy group 
of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group pattern (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"

    # Pattern for decanoyl group (10-carbon chain attached to ester)
    # More flexible pattern allowing for substitutions
    decanoyl_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
    
    if not mol.HasSubstructMatch(decanoyl_pattern):
        return False, "No 10-carbon chain attached to ester group found"
    
    # Get matches
    matches = mol.GetSubstructMatches(decanoyl_pattern)
    
    for match in matches:
        # Get the carbon atoms in the chain
        chain_atoms = match[2:]  # Skip O-C(=O)
        
        # Check if chain is linear (each carbon should have 2 or 3 non-H neighbors except terminal)
        is_valid_chain = True
        for i, atom_idx in enumerate(chain_atoms):
            atom = mol.GetAtomWithIdx(atom_idx)
            non_h_neighbors = len([n for n in atom.GetNeighbors() 
                                 if n.GetAtomicNum() != 1 and n.GetIdx() not in chain_atoms])
            
            # Terminal carbon can have 1 non-chain neighbor
            if i == len(chain_atoms)-1:
                if non_h_neighbors > 1:
                    is_valid_chain = False
                    break
            # Other carbons should have at most 1 non-chain neighbor
            elif non_h_neighbors > 1:
                is_valid_chain = False
                break
        
        if is_valid_chain:
            # Verify the ester oxygen is connected to carbon (not part of phosphate, etc)
            o_atom = mol.GetAtomWithIdx(match[0])
            neighbors = [x for x in o_atom.GetNeighbors() if x.GetIdx() != match[1]]
            if len(neighbors) != 1:
                continue
            if neighbors[0].GetAtomicNum() != 6:
                continue
                
            # Molecular weight sanity check (should be > 170 Da for smallest decanoate ester)
            mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
            if mol_wt < 170:
                return False, "Molecular weight too low for decanoate ester"
                
            return True, "Contains decanoyl group connected via ester linkage with 10-carbon chain"
    
    return False, "No valid decanoate ester group found"