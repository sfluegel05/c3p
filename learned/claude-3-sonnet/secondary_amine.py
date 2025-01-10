"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:32854 secondary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains a secondary amine group.
    A secondary amine has exactly one N-H bond and two non-hydrogen substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Patterns to exclude
    exclude_patterns = {
        "amide": "[NX3H1]C(=O)",
        "sulfonamide": "[NX3H1]S(=O)(=O)",
        "nitro": "[NX3]([O-])=O",
        "nitroso": "[NX2]=O",
        "imine": "[NX2]=C",
        "azide": "[N-][N+]#N",
        "aromatic_n": "n",
        "guanidine": "[NX3H1]C(=N)N",
        "urea": "[NX3H1]C(=O)N",
        "carbamate": "[NX3H1]C(=O)O",
        "charged_n": "[$([N+]),$([N-])]",
        "thiourea": "[NX3H1]C(=S)N",
        "phosphoramide": "[NX3H1]P(=O)",
        "n_oxide": "[NX4+]",
    }

    # Convert patterns to RDKit molecules
    exclude_mols = {name: Chem.MolFromSmarts(pattern) 
                   for name, pattern in exclude_patterns.items()}

    # Find all nitrogens with exactly one hydrogen
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3H1]"))
    
    for match in matches:
        N_idx = match[0]
        N_atom = mol.GetAtomWithIdx(N_idx)
        
        # Skip if nitrogen is aromatic
        if N_atom.GetIsAromatic():
            continue
            
        # Count non-hydrogen neighbors
        heavy_neighbors = len([n for n in N_atom.GetNeighbors() 
                             if n.GetAtomicNum() != 1])
        
        if heavy_neighbors != 2:
            continue
            
        # Check if this nitrogen is part of any excluded patterns
        excluded = False
        for name, pattern in exclude_mols.items():
            if pattern is not None and mol.HasSubstructMatch(pattern):
                matches = mol.GetSubstructMatches(pattern)
                for m in matches:
                    if N_idx in m:
                        excluded = True
                        break
            if excluded:
                break
                
        if not excluded:
            # Additional checks for specific cases
            neighbors = [n for n in N_atom.GetNeighbors() if n.GetAtomicNum() != 1]
            
            # Check if both neighbors are carbons (most common case)
            if all(n.GetAtomicNum() == 6 for n in neighbors):
                return True, "Contains secondary amine group (NH with two carbon substituents)"
            
            # Check if neighbors are C, O, or alkyl groups
            valid_neighbors = True
            for n in neighbors:
                if n.GetAtomicNum() not in [6, 8]:  # C or O
                    valid_neighbors = False
                elif n.GetIsAromatic():
                    # Allow aromatic carbons
                    continue
                elif any(bond.GetBondType() != Chem.BondType.SINGLE 
                        for bond in n.GetBonds()):
                    valid_neighbors = False
                    
            if valid_neighbors:
                return True, "Contains secondary amine group (NH with two non-H substituents)"

    return False, "No valid secondary amine group found"