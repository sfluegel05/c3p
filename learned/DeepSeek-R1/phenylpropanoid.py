"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids have a phenylpropane skeleton (C6-C3) with possible substitutions.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phenylpropane core: benzene connected to three carbons (C6-C3)
    # SMARTS allows for any bonds and substituents on the three carbons
    core_pattern = Chem.MolFromSmarts('[c]1[c][c][c][c][c]1-[C]~[C]~[C]')
    if not mol.HasSubstructMatch(core_pattern):
        # Alternative check with possible oxygen bridge (e.g., esters)
        # Looks for benzene-O-C-C-C or benzene-C-C-C with substitutions
        alt_core_pattern = Chem.MolFromSmarts('[c]1[c][c][c][c][c]1-[O]-[C]~[C]~[C]')
        if not mol.HasSubstructMatch(alt_core_pattern):
            return False, "No phenylpropane skeleton (C6-C3) detected"

    # Verify presence of aromatic system (at least one aromatic ring)
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) 
                      if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if not aromatic_rings:
        return False, "No aromatic rings found"

    # Optional: Check for common oxygen-containing groups (but not strictly required)
    oxygen_pattern = Chem.MolFromSmarts('[O]')
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No oxygen atoms found (common but not mandatory)"

    return True, "Contains phenylpropane skeleton with aromatic features"