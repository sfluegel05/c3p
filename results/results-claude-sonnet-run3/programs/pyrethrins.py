from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrethrins(smiles: str):
    """
    Determines if a molecule belongs to the pyrethrin class of insecticides.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a pyrethrin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Core structure with cyclopentenone ring and ester linkage
    core_pattern = Chem.MolFromSmarts('[#6]-1-[#6]-[#6](=[O])-[#6](-[#6])=[#6]-1-[#6]')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing cyclopentenone core"

    # Ester linkage
    ester = Chem.MolFromSmarts('[#6]-[#6](=[O])-[#8]-[#6]')
    if not mol.HasSubstructMatch(ester):
        return False, "Missing ester linkage"

    # Cyclopropane ring
    cyclopropane = Chem.MolFromSmarts('[#6]-1-[#6]-[#6]-1')
    if not mol.HasSubstructMatch(cyclopropane):
        return False, "Missing cyclopropane ring"

    # Specific pattern for pyrethrins (more detailed core structure)
    pyrethrin_core = Chem.MolFromSmarts('CC1=C(C)C(OC(=O)C2C(C=C(C)C)C2(C)C)CC1=O')
    if not mol.HasSubstructMatch(pyrethrin_core):
        return False, "Missing characteristic pyrethrin core structure"

    # Check for characteristic side chains
    type_2_ester = Chem.MolFromSmarts('COC(=O)C(C)=C')
    
    if mol.HasSubstructMatch(type_2_ester):
        if "C\\C=C/C=C" in smiles:
            return True, "Pyrethrin II"
        elif "C\\C=C/C" in smiles:
            return True, "Cinerin II"
        elif "CC\\C=C/" in smiles:
            return True, "Jasmolin II"
    else:
        if "C\\C=C/C=C" in smiles:
            return True, "Pyrethrin I"
        elif "C\\C=C/C" in smiles:
            return True, "Cinerin I"
        elif "CC\\C=C/" in smiles:
            return True, "Jasmolin I"

    return True, "Unspecified pyrethrin"
# Pr=1.0
# Recall=1.0