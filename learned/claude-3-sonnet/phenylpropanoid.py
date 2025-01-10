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
        # Basic phenylpropane skeleton (more flexible)
        "phenylpropane": "[$(c1ccccc1CCC),$(c1ccccc1CC=C),$(c1ccccc1C(=O)CC),$(c1ccccc1CCO)]",
        
        # Coumarin core (more general)
        "coumarin": "[$(O=C1Oc2ccccc2C=C1),$(O=C1Oc2ccccc2CC1)]",
        
        # Chromone core
        "chromone": "[$(O=C1C=COc2ccccc12),$(O=C1CCOc2ccccc12)]",
        
        # Benzofuran core
        "benzofuran": "c1ccc2occc2c1",
        
        # Flavonoid cores (more comprehensive)
        "flavonoid": "[$(O=C1CC(c2ccccc2)Oc2ccccc12),$(O=C1C(c2ccccc2)=COc2ccccc12)]",
        
        # Isoflavonoid core (more flexible)
        "isoflavonoid": "[$(O=C1C=C(c2ccccc2)Oc2ccccc12),$(O1C=C(c2ccccc2)C(=O)c2ccccc12)]",
        
        # Lignin building blocks (monolignols)
        "monolignol": "[$(c1cc(O)c(O)cc1CC=CO),$(c1cc(OC)c(O)cc1CC=CO),$(c1cc(O)cc(OC)c1CC=CO)]",
        
        # Cinnamate derivatives
        "cinnamate": "[$(C=CC(=O)O),$(C=CC(=O)OC),$(C=CC(=O)N)]",
        
        # General oxygen-containing heterocycle pattern
        "oxyheterocycle": "[$(O=C1Oc2ccccc2C=C1),$(O=C1Oc2ccccc2CC1),$(c1cc2OCOc2cc1)]"
    }

    # Check for any of the characteristic substructures
    found_substructures = []
    for name, smarts in substructures.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_substructures.append(name)

    # Additional structure checks
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for x in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in x))
    
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Decision logic
    if found_substructures:
        reason = f"Contains phenylpropanoid features: {', '.join(found_substructures)}. "
        reason += f"Has {aromatic_rings} aromatic rings and {o_count} oxygen atoms."
        return True, reason
    
    # Check for complex polycyclic structures with oxygen
    if aromatic_rings >= 2 and o_count >= 2 and mol_wt >= 250:
        # Look for characteristic oxygen patterns in aromatic context
        o_pattern = Chem.MolFromSmarts("[$(cOc),$(C(=O)Oc),$(COc)]")
        if mol.HasSubstructMatch(o_pattern):
            return True, f"Complex polycyclic structure with {aromatic_rings} aromatic rings and {o_count} oxygens in characteristic patterns"

    return False, "No characteristic phenylpropanoid substructures found"