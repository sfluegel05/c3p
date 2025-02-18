"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: CHEBI:28700 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide (cerebroside with galactose).
    Structure must have:
    - Galactose (pyranose form, any substitution pattern)
    - Connected via glycosidic bond to ceramide:
      - Sphingoid base (long chain with amino and hydroxyl groups)
      - Amide-linked fatty acid chain (long aliphatic chain)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Match galactose core (pyranose ring with at least 4 oxygen substituents, allowing modifications)
    # Pattern accounts for possible substitutions (e.g. sulfates) on hydroxyl groups
    galactose_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@@H]([C@H]([C@@H](O1)CO[*:1])[*:2])[*:3])[*:4]")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose moiety detected"

    # Find glycosidic bond between galactose and ceramide
    # Look for O connecting galactose's anomeric carbon to another atom
    anomeric_carbon = Chem.MolFromSmarts("[C@H](-[O]-[#6])-[C@H]")  # Anomeric C-O-C link
    if not mol.HasSubstructMatch(anomeric_carbon):
        return False, "No glycosidic bond detected"

    # Check for ceramide structure: sphingosine + fatty acid amide
    # Sphingosine base pattern: long chain with adjacent NH and OH groups
    sphingo_pattern = Chem.MolFromSmarts("[C]([NH][C](=O)[C][C])[C][C][C@H](O)[C@H]")  # Approx sphingosine backbone
    if not mol.HasSubstructMatch(sphingo_pattern):
        return False, "No sphingosine base detected"

    # Verify amide-linked fatty acid (long chain, at least 14 carbons)
    amide_n = Chem.MolFromSmarts("[NX3][C](=O)[CX4][CX4,H2]")
    amide_matches = mol.GetSubstructMatches(amide_n)
    if not amide_matches:
        return False, "No amide bond found"

    # Check fatty acid chain length from amide
    for match in amide_matches:
        fatty_start = mol.GetAtomWithIdx(match[2])  # Carbon adjacent to amide
        chain = AllChem.FindAtomEnvironmentOfRadiusN(mol, 20, fatty_start.GetIdx())
        chain_atoms = {mol.GetBondWithIdx(i).GetBeginAtomIdx() for i in chain} | \
                     {mol.GetBondWithIdx(i).GetEndAtomIdx() for i in chain}
        carbon_count = sum(1 for idx in chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if carbon_count >= 14:
            return True, "Galactose linked to ceramide with long-chain fatty acid"

    return False, "Insufficient fatty acid chain length"