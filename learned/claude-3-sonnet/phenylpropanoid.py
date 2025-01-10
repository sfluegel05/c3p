"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26171 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids are compounds with a phenyl ring and a 3-carbon propane skeleton,
    including their derivatives like flavonoids, coumarins, and lignins.

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

    # Look for aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No phenyl (aromatic) ring found"

    # Common phenylpropanoid substructures
    substructures = {
        # Basic phenylpropane skeleton
        "phenylpropane": "[$(c1ccccc1CCC),$(c1ccccc1CC=C),$(c1ccccc1C(=O)CC)]",
        
        # Coumarin core
        "coumarin": "O=C1OC=Cc2ccccc12",
        
        # Flavonoid core (includes flavones, flavanones, etc.)
        "flavonoid": "O=C1CC(c2ccccc2)Oc2ccccc12",
        
        # Isoflavonoid core
        "isoflavonoid": "O=C1C=C(c2ccccc2)Oc2ccccc12",
        
        # Lignin building blocks (monolignols)
        "monolignol": "[$(c1cc(O)c(O)cc1CC=CO),$(c1cc(OC)c(O)cc1CC=CO)]",
        
        # Cinnamate derivatives
        "cinnamate": "[$(C=CC(=O)O),$(C=CC(=O)OC)]",
    }

    # Check for any of the characteristic substructures
    found_substructures = []
    for name, smarts in substructures.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_substructures.append(name)

    if not found_substructures:
        return False, "No characteristic phenylpropanoid substructures found"

    # Additional checks
    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for x in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in x))
    
    # Count oxygen atoms (phenylpropanoids typically have oxygen-containing groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Most phenylpropanoids have at least one oxygen
    if o_count == 0:
        return False, "No oxygen atoms found - unusual for phenylpropanoids"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Very small molecules are unlikely to be phenylpropanoids
    if mol_wt < 120:  # Approximate minimum weight for phenylpropanoids
        return False, "Molecule too small to be a phenylpropanoid"

    reason = f"Contains phenylpropanoid features: {', '.join(found_substructures)}. "
    reason += f"Has {aromatic_rings} aromatic rings and {o_count} oxygen atoms."
    
    return True, reason