"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids are organic aromatic compounds with a structure based on a
    phenylpropane skeleton (C6-C3 unit). This includes subclasses such as flavonoids,
    anthocyanins, coumarins, lignans, and stilbenes.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for phenylpropanoid subclasses

    # General phenylpropane skeleton: benzene ring attached to a C3 chain
    phenylpropane_pattern = Chem.MolFromSmarts("c1ccccc1CCC")
    phenylpropenoid_pattern = Chem.MolFromSmarts("c1ccccc1C=CC")  # Includes double bonds
    cinnamic_acid_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)O")
    cinnamaldehyde_pattern = Chem.MolFromSmarts("c1ccccc1C=CC=O")
    cinnamyl_alcohol_pattern = Chem.MolFromSmarts("c1ccccc1C=CCO")

    # Flavonoids: 15-carbon skeleton with two benzene rings and a heterocyclic ring
    flavonoid_pattern = Chem.MolFromSmarts("c1cc(-c2ccc3occ(c3c2)c2ccccc2)ccc1")
    
    # Coumarins: benzopyran-2-one core structure
    coumarin_pattern = Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2O1")
    
    # Lignans: dimers derived from phenylpropanoids, C6-C3 units linked
    lignan_pattern = Chem.MolFromSmarts("c1ccccc1C[C@@H](C)c2cccc(O)c2")
    
    # Anthocyanins: flavylium ion core
    anthocyanin_pattern = Chem.MolFromSmarts("[O+]c1cc2ccccc2oc1")
    
    # Stilbenes: 1,2-diphenylethylene core
    stilbene_pattern = Chem.MolFromSmarts("c1ccc(cc1)/C=C/c2ccccc2")
    
    # Chalcones: open-chain flavonoids with two aromatic rings linked by a three-carbon α,β-unsaturated carbonyl system
    chalcone_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)c1ccccc1")
    
    # Phenylpropanoid glycosides: phenylpropanoid attached to sugar moiety
    glycoside_pattern = Chem.MolFromSmarts("c1ccccc1C=CCO[C;X4]")  # O-glycosides
    
    # List of patterns with descriptions
    patterns = [
        ("phenylpropane skeleton", phenylpropane_pattern),
        ("phenylpropenoid skeleton", phenylpropenoid_pattern),
        ("cinnamic acid structure", cinnamic_acid_pattern),
        ("cinnamaldehyde structure", cinnamaldehyde_pattern),
        ("cinnamyl alcohol structure", cinnamyl_alcohol_pattern),
        ("flavonoid core structure", flavonoid_pattern),
        ("coumarin core structure", coumarin_pattern),
        ("lignan structure", lignan_pattern),
        ("anthocyanin core structure", anthocyanin_pattern),
        ("stilbene core structure", stilbene_pattern),
        ("chalcone core structure", chalcone_pattern),
        ("phenylpropanoid glycoside", glycoside_pattern),
    ]

    # Check for matches with specific patterns
    for description, pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {description}"

    # Additional checks for oxygenated functionalities common in phenylpropanoids
    phenylpropanoid_like = False
    ring_info = mol.GetRingInfo()
    num_aromatic_rings = sum(1 for ring in ring_info.BondRings() if all(mol.GetBondWithIdx(idx).GetIsAromatic() for idx in ring))
    
    # Check if molecule has at least one aromatic ring
    if num_aromatic_rings >= 1:
        # Look for C6-C3 unit (total carbons >=9)
        num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if num_carbons >= 9:
            # Look for hydroxyl or methoxy groups attached to aromatic ring
            phenolic_oxygen = Chem.MolFromSmarts("c[OH]")
            methoxy_group = Chem.MolFromSmarts("cOC")
            if mol.HasSubstructMatch(phenolic_oxygen) or mol.HasSubstructMatch(methoxy_group):
                phenylpropanoid_like = True
                
    if phenylpropanoid_like:
        return True, "Contains phenolic aromatic ring with aliphatic chain"

    return False, "Does not match phenylpropanoid structures"