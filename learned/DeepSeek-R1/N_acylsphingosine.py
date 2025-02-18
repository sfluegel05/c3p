"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:17578 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    Must contain:
    1. Sphingosine backbone with amino alcohol (NH connected to carbonyl and adjacent OH group)
    2. At least 12 carbons in the sphingosine chain
    3. At least one double bond in the sphingosine chain
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find amide-linked amino alcohol pattern: N-C(=O)-C(OH)-...
    # SMARTS: [NX3][CX3](=[OX1])[CX4H1][OX2H1]
    amide_amino_alcohol = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4H1][OX2H1]")
    matches = mol.GetSubstructMatches(amide_amino_alcohol)
    if not matches:
        return False, "No amide-linked amino alcohol found"
    
    for match in matches:
        nitrogen_idx = match[0]
        sphingosine_carbon_idx = match[2]
        sphingosine_carbon = mol.GetAtomWithIdx(sphingosine_carbon_idx)
        
        # Traverse sphingosine chain (main chain only)
        chain_length = 1  # include current carbon
        double_bond_found = False
        current = sphingosine_carbon
        prev_atom = mol.GetAtomWithIdx(match[1])  # carbonyl carbon (don't backtrack)
        
        while True:
            next_carbons = []
            for neighbor in current.GetNeighbors():
                if neighbor.GetAtomicNum() != 6:
                    continue
                if neighbor.GetIdx() == prev_atom.GetIdx():
                    continue
                next_carbons.append(neighbor)
            
            # Select next carbon in main chain (prioritize non-branching)
            if not next_carbons:
                break
            next_atom = max(next_carbons, key=lambda a: a.GetDegree(), default=None)
            if not next_atom:
                break
            
            # Check bond for double bond
            bond = current.GetBondBetweenAtoms(next_atom.GetIdx(), current.GetIdx())
            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                double_bond_found = True
            
            chain_length += 1
            prev_atom = current
            current = next_atom
        
        # Check chain requirements
        if chain_length >= 12 and double_bond_found:
            return True, "Contains sphingosine backbone with N-linked fatty acyl chain"
    
    return False, "Missing sphingosine backbone features or chain too short"