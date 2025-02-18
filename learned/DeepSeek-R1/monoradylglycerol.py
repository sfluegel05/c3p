"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: CHEBI monoradylglycerol (Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone with exactly one acyl (ester), alkyl (ether), or alk-1-enyl (vinyl ether) substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define glycerol backbone patterns with substituent in any position
    glycerol_patterns = [
        Chem.MolFromSmarts("[CH2](-O[!H0])-[CH](-O)-[CH2](-O)"),  # substituent on C1
        Chem.MolFromSmarts("[CH2](-O)-[CH](-O[!H0])-[CH2](-O)"),  # substituent on C2
        Chem.MolFromSmarts("[CH2](-O)-[CH](-O)-[CH2](-O[!H0])")   # substituent on C3
    ]

    # Check if any glycerol pattern matches
    backbone_found = any(mol.HasSubstructMatch(patt) for patt in glycerol_patterns if patt)
    if not backbone_found:
        return False, "No glycerol backbone with two hydroxyls and one O-linked group"

    # Find all oxygen atoms in non-hydroxyl linkages (substituents)
    substituent_oxygens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and not atom.GetTotalNumHs():
            # Check if connected to glycerol carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    # Verify neighbor is part of glycerol backbone
                    for patt in glycerol_patterns:
                        if patt and neighbor.HasSubstructMatch(patt):
                            substituent_oxygens.append(atom)
                            break

    if len(substituent_oxygens) != 1:
        return False, f"Found {len(substituent_oxygens)} substituents, need exactly 1"

    # Check substituent type using SMARTS patterns
    ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2][CX3]=[OX1]"))  # Acyl (ester)
    ether = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2][CX4H2]"))      # Alkyl (ether)
    vinyl_ether = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2]C=C"))    # Alk-1-enyl (vinyl ether)

    if not any([ester, ether, vinyl_ether]):
        return False, "Substituent is not acyl, alkyl, or alk-1-enyl type"

    # Determine exact substituent type
    if ester:
        return True, "Glycerol with one acyl (ester) substituent"
    elif vinyl_ether:
        return True, "Glycerol with one alk-1-enyl (vinyl ether) substituent"
    else:
        return True, "Glycerol with one alkyl (ether) substituent"