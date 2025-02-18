"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: glycosphingolipid (CHEBI:24404)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    Criteria: Carbohydrate attached via glycosidic linkage to ceramide/sphingoid base.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Combined pattern: sphingoid/ceramide core connected via glycosidic O to a sugar
    # 1. Ceramide part: amide group (N-C=O) attached to long chain
    # 2. Glycosidic linkage: O connecting to a sugar (ring with multiple OH/O)
    # 3. Sugar: at least 4 contiguous C-O/C-OH in a ring (heuristic for carbohydrate)
    pattern = Chem.MolFromSmarts(
        "[NX3][C](=[OX1])[CX4]."  # Ceramide amide group
        "[OX2][C@H]("  # Glycosidic oxygen connected to sphingoid carbon
            "[C@H]([OH0])[C@H](O)[CH2,CH][OH0]"  # Sphingoid backbone (simplified)
        ")-!@*"  # Disconnected from sugar part (avoid same molecule)
        ".[O;R]-[C;R]-[C;R]-[C;R]-[C;R]-[O;R]|"  # Pyranose ring (O in ring)
        "[O;R]-[C;R]1-[C;R]-[C;R]-[C;R]-[O;R]-1"  # Furanose ring
    )
    
    if not mol.HasSubstructMatch(pattern):
        # Fallback: Check for any O-linked sugar to ceramide/sphingoid
        # General glycosphingolipid pattern: ceramide/sphingoid + O-linked sugar
        general_pattern = Chem.MolFromSmarts(
            "[NX3,CX3]([CX4])!@[OX2][C]("  # O linking to NH/CO and sugar
            ").([#6]1@[#6](-[OH])@[#6](-[OH])@[#6]@[#6]@[#6]@1)"  # Carbohydrate ring
        )
        if not mol.HasSubstructMatch(general_pattern):
            return False, "Missing glycosidic linkage between sugar and ceramide/sphingoid"
    
    return True, "Carbohydrate attached via glycosidic bond to ceramide/sphingoid"