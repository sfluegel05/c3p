"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: CHEBI:38836 anthoxanthins
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid pigments in plants, ranging in color from white/colorless
    to yellow. They are water-soluble and often found in flower petals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for common anthoxanthin backbones
    backbones = ['c1c(o[c@@H]2[c@H](c3ccc(O)cc3o2)ccc1O',  # flavonol (quercetin)
                 'c1c(oc2cc(O)cc(O)c2c1=O)',  # flavone
                 'C=C(C=Cc1ccc(O)c(O)c1)C(=O)/C=C/c1ccc(O)cc1',  # chalcone
                 'c1c(oc2ccc(O)cc2c1=O)',  # flavanone
                 'c1c(-c2ccc(O)cc2)oc2cc(O)cc(O)c2c1=O',  # isoflavone
                 # Add more patterns as needed
                 ]
    has_backbone = any(mol.HasSubstructMatch(Chem.MolFromSmarts(b)) for b in backbones)
    if not has_backbone:
        return False, "No common anthoxanthin backbone found"

    # Check for common anthoxanthin substituents
    substituents = ['[OX2H]',  # hydroxyl
                    '[OX2][CX4][CX4]',  # glycosidic linkage
                    '[CX3](=O)[OX2]',  # ester/carboxylate
                    '[SX4+2]',  # sulfate
                    '[CX3][OX2][CH3]',  # methoxy
                    '[CX3](-[OX2][CX3]=O)-[CX3]=O',  # malonyl
                    # Add more patterns as needed
                    ]
    has_substituents = any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in substituents)
    if not has_substituents:
        return False, "No common anthoxanthin substituents found"

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    hba = rdMolDescriptors.CalcNumLipinskiHBA(mol)
    hbd = rdMolDescriptors.CalcNumLipinskiHBD(mol)

    if mol_wt < 150 or mol_wt > 1500:  # Adjusted ranges
        return False, "Molecular weight out of typical anthoxanthin range"
    if logp > 7:  # Adjusted range
        return False, "LogP too high for water-soluble anthoxanthin"
    if hba < 3 or hbd < 1:  # Adjusted ranges
        return False, "Insufficient hydrogen bonding for water solubility"

    return True, "Contains anthoxanthin backbone and substituents, with properties consistent with water-soluble pigments"