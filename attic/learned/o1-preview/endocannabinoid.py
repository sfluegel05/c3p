"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.

    Endocannabinoids are endogenous lipid-based molecules that bind to cannabinoid receptors.
    They typically have:
    - A long-chain fatty acid chain (often with multiple double bonds)
    - A head group like ethanolamine or glycerol
    - Linked via an amide or ester bond

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify amide and ester bonds
    amide_bond_smarts = "[NX3][CX3](=O)[C]"
    ester_bond_smarts = "[OX2H0][CX3](=O)[C]"
    amide_bond = Chem.MolFromSmarts(amide_bond_smarts)
    ester_bond = Chem.MolFromSmarts(ester_bond_smarts)

    # Find amide and ester bonds
    amide_matches = mol.GetSubstructMatches(amide_bond)
    ester_matches = mol.GetSubstructMatches(ester_bond)

    # Collect potential linkages
    linkages = []
    for match in amide_matches:
        linkages.append(('amide', match))
    for match in ester_matches:
        linkages.append(('ester', match))

    if not linkages:
        return False, "No amide or ester linkage found"

    # Define head group patterns
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")

    # For each linkage, check for head group and fatty acid chain
    for linkage_type, match in linkages:
        if linkage_type == 'amide':
            nitrogen_idx = match[0]
            carbon_idx = match[2]
            bond_idx = mol.GetBondBetweenAtoms(nitrogen_idx, match[1]).GetIdx()
            # Fragment molecule at the amide bond
            fragments = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True), asMols=True)
            # Assign fragments to head group and fatty acid chain
            headgroup_frag = None
            fattyacid_frag = None
            for frag in fragments:
                if frag.HasSubstructMatch(ethanolamine_pattern):
                    headgroup_frag = frag
                else:
                    fattyacid_frag = frag
            if headgroup_frag and fattyacid_frag:
                # Check length of fatty acid chain
                num_carbons = sum(1 for atom in fattyacid_frag.GetAtoms() if atom.GetAtomicNum() == 6)
                if num_carbons >= 12:
                    return True, "Endocannabinoid with ethanolamine head group and fatty acid chain found"
            else:
                continue  # Try next linkage

        elif linkage_type == 'ester':
            oxygen_idx = match[0]
            carbon_idx = match[2]
            bond_idx = mol.GetBondBetweenAtoms(oxygen_idx, match[1]).GetIdx()
            # Fragment molecule at the ester bond
            fragments = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True), asMols=True)
            # Assign fragments to head group and fatty acid chain
            headgroup_frag = None
            fattyacid_frag = None
            for frag in fragments:
                if frag.HasSubstructMatch(glycerol_pattern):
                    headgroup_frag = frag
                else:
                    fattyacid_frag = frag
            if headgroup_frag and fattyacid_frag:
                # Check length of fatty acid chain
                num_carbons = sum(1 for atom in fattyacid_frag.GetAtoms() if atom.GetAtomicNum() == 6)
                if num_carbons >= 12:
                    return True, "Endocannabinoid with glycerol head group and fatty acid chain found"
            else:
                continue  # Try next linkage

    return False, "Molecule does not match endocannabinoid criteria"