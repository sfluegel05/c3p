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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern with atom indices
    glycerol_pattern = Chem.MolFromSmarts('[C:1]-[C:2]-[C:3]')

    # Find glycerol backbone matches
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    found_lysophosphatidic_acid = False

    for match in matches:
        c1_idx, c2_idx, c3_idx = match
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Check for hydroxyl group on central carbon (c2)
        has_hydroxyl_c2 = False
        for neighbor in c2.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                # Hydroxyl oxygen
                has_hydroxyl_c2 = True
                break
        if not has_hydroxyl_c2:
            continue  # Move to next match if no hydroxyl on c2

        # Option 1: c1 connected to phosphate, c3 connected to acyl chain
        is_c1_phosphate = False
        for bond in c1.GetBonds():
            nbr = bond.GetOtherAtom(c1)
            if nbr.GetAtomicNum() == 8:
                oxygen = nbr
                for obond in oxygen.GetBonds():
                    onbr = obond.GetOtherAtom(oxygen)
                    if onbr.GetAtomicNum() == 15:
                        is_c1_phosphate = True
                        break
                if is_c1_phosphate:
                    break

        is_c3_acyl_chain = False
        for bond in c3.GetBonds():
            nbr = bond.GetOtherAtom(c3)
            if nbr.GetAtomicNum() == 8:
                oxygen = nbr
                # Check for ester linkage to acyl chain
                for obond in oxygen.GetBonds():
                    onbr = obond.GetOtherAtom(oxygen)
                    if onbr.GetAtomicNum() == 6 and onbr.GetIdx() != c3_idx:
                        carbonyl_carbon = onbr
                        # Check for carbonyl group (C=O)
                        is_carbonyl = False
                        for cbond in carbonyl_carbon.GetBonds():
                            cnbr = cbond.GetOtherAtom(carbonyl_carbon)
                            if cnbr.GetAtomicNum() == 8 and cbond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                is_carbonyl = True
                                break
                        if is_carbonyl:
                            is_c3_acyl_chain = True
                            break
                if is_c3_acyl_chain:
                    break

        if is_c1_phosphate and is_c3_acyl_chain:
            found_lysophosphatidic_acid = True
            break

        # Option 2: c3 connected to phosphate, c1 connected to acyl chain
        is_c3_phosphate = False
        for bond in c3.GetBonds():
            nbr = bond.GetOtherAtom(c3)
            if nbr.GetAtomicNum() == 8:
                oxygen = nbr
                for obond in oxygen.GetBonds():
                    onbr = obond.GetOtherAtom(oxygen)
                    if onbr.GetAtomicNum() == 15:
                        is_c3_phosphate = True
                        break
                if is_c3_phosphate:
                    break

        is_c1_acyl_chain = False
        for bond in c1.GetBonds():
            nbr = bond.GetOtherAtom(c1)
            if nbr.GetAtomicNum() == 8:
                oxygen = nbr
                # Check for ester linkage to acyl chain
                for obond in oxygen.GetBonds():
                    onbr = obond.GetOtherAtom(oxygen)
                    if onbr.GetAtomicNum() == 6 and onbr.GetIdx() != c1_idx:
                        carbonyl_carbon = onbr
                        # Check for carbonyl group (C=O)
                        is_carbonyl = False
                        for cbond in carbonyl_carbon.GetBonds():
                            cnbr = cbond.GetOtherAtom(carbonyl_carbon)
                            if cnbr.GetAtomicNum() == 8 and cbond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                is_carbonyl = True
                                break
                        if is_carbonyl:
                            is_c1_acyl_chain = True
                            break
                if is_c1_acyl_chain:
                    break

        if is_c3_phosphate and is_c1_acyl_chain:
            found_lysophosphatidic_acid = True
            break

    if found_lysophosphatidic_acid:
        return True, "Molecule is a lysophosphatidic acid"
    else:
        return False, "Does not match lysophosphatidic acid structure"