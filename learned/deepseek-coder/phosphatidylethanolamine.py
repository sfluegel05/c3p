"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: CHEBI:17544 phosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is a glycerophospholipid with a phosphatidyl group esterified to ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached to fatty acids and 1 oxygen attached to phosphate)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for phosphate group attached to glycerol
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"

    # Check if phosphate is attached to glycerol
    phosphate_attached = False
    for match in phosphate_matches:
        phosphate_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in phosphate_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and mol.GetSubstructMatch(glycerol_pattern):
                phosphate_attached = True
                break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "Phosphate group not attached to glycerol backbone"

    # More flexible ethanolamine pattern to catch different protonation states and connectivity
    ethanolamine_pattern = Chem.MolFromSmarts("[NX3,nX3+][CH2X4][OX2,oX1-]?")
    ethanolamine_matches = mol.GetSubstructMatches(ethanolamine_pattern)
    if len(ethanolamine_matches) == 0:
        return False, "No ethanolamine group found"

    # Check if ethanolamine is attached to phosphate (directly or through oxygen)
    ethanolamine_attached = False
    for match in ethanolamine_matches:
        ethanolamine_atom = mol.GetAtomWithIdx(match[0])
        # Check direct P-N bond or P-O-C-N connection
        for neighbor in ethanolamine_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'P':
                ethanolamine_attached = True
                break
            elif neighbor.GetSymbol() == 'O':
                for o_neighbor in neighbor.GetNeighbors():
                    if o_neighbor.GetSymbol() == 'P':
                        ethanolamine_attached = True
                        break
                if ethanolamine_attached:
                    break
        if ethanolamine_attached:
            break
    if not ethanolamine_attached:
        return False, "Ethanolamine group not attached to phosphate"

    # Look for two ester groups (-O-C(=O)-) attached to glycerol
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - phosphatidylethanolamines typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylethanolamine"

    return True, "Contains glycerol backbone with two fatty acid chains, a phosphate group, and an ethanolamine group"