"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a monoacylglycerol phosphate, consisting of a glycerol backbone
    esterified with one fatty acid chain and phosphorylated on one of the hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import AllChem

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find glycerol backbone (three connected carbons)
    # Define SMARTS pattern for glycerol backbone carbons
    glycerol_pattern = Chem.MolFromSmarts('C-C-C')
    matches = mol.GetSubstructMatches(glycerol_pattern)

    if not matches:
        return False, "No glycerol backbone found"

    found_lysophosphatidic_acid = False
    for match in matches:
        c1_idx, c2_idx, c3_idx = match

        # Now, check the substituents of each carbon
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Initialize substituent flags
        has_phosphate = False
        has_ester = False
        has_hydroxyl = False

        # Check substituents of c1 for phosphate group
        for bond in c1.GetBonds():
            nbr = bond.GetOtherAtom(c1)
            if nbr.GetAtomicNum() == 8:  # Oxygen atom
                # Check if this oxygen is connected to a phosphate group
                oxygen = nbr
                for obond in oxygen.GetBonds():
                    phosphorus = obond.GetOtherAtom(oxygen)
                    if phosphorus.GetAtomicNum() == 15:  # Phosphorus atom
                        has_phosphate = True

        # Check substituents of c2 for hydroxyl group
        for bond in c2.GetBonds():
            nbr = bond.GetOtherAtom(c2)
            if nbr.GetAtomicNum() == 8:  # Oxygen atom
                if nbr.GetDegree() == 1:
                    has_hydroxyl = True

        # Check substituents of c3 for ester group
        for bond in c3.GetBonds():
            nbr = bond.GetOtherAtom(c3)
            if nbr.GetAtomicNum() == 8:  # Oxygen atom
                oxygen = nbr
                for obond in oxygen.GetBonds():
                    carbonyl = obond.GetOtherAtom(oxygen)
                    if carbonyl.GetAtomicNum() == 6 and carbonyl.GetIdx() != c3_idx:
                        # Check if carbon is a carbonyl carbon (C=O)
                        is_carbonyl = False
                        for cbond in carbonyl.GetBonds():
                            o_nbr = cbond.GetOtherAtom(carbonyl)
                            if o_nbr.GetAtomicNum() == 8 and o_nbr.GetIdx() != oxygen.GetIdx():
                                is_carbonyl = True
                        if is_carbonyl:
                            has_ester = True

        # Check if we have all required substituents
        if has_phosphate and has_hydroxyl and has_ester:
            found_lysophosphatidic_acid = True
            break

    if found_lysophosphatidic_acid:
        return True, "Molecule is a lysophosphatidic acid"
    else:
        return False, "Does not match lysophosphatidic acid structure"