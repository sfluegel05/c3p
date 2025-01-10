"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is defined as 4-(2-aminoethyl)pyrocatechol and derivatives formed by substitution.

    Catecholamines have a catechol moiety (benzene ring with adjacent hydroxyls)
    and an aminoethyl side chain attached to the ring at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the catecholamine core structure
    # This pattern matches a benzene ring with hydroxyl groups on adjacent carbons (positions 1 and 2)
    # and an aminoethyl chain attached to the ring at position 4 relative to the hydroxyl at position 1.
    catecholamine_pattern = Chem.MolFromSmarts("""
    [$([OH]),$(O-[H])]c1ccc([CX2][CX3][NX3H2,NX3H,NX3+])c([OH])c1
    """)

    if mol.HasSubstructMatch(catecholamine_pattern):
        return True, "Molecule matches the catecholamine core structure"

    # If not matched, provide detailed reasoning
    # Check for catechol moiety
    catechol_pattern = Chem.MolFromSmarts('c1c(O)ccc(O)c1')
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "Molecule does not contain a catechol moiety (adjacent dihydroxybenzene)"

    # Check for aminoethyl side chain attached to the ring
    # Define aminoethyl group pattern (allowing for substitutions on nitrogen)
    aminoethyl_pattern = Chem.MolFromSmarts('[CX3][CX3][NX3]')
    # Find attachment points of the aminoethyl group
    attachment_found = False
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        # Check if bond connects the ring to an aminoethyl group
        if begin_atom.IsInRing() != end_atom.IsInRing():
            # Get the non-ring atom
            side_chain_atom = begin_atom if not begin_atom.IsInRing() else end_atom
            # Get the environment around the side chain atom
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, 2, side_chain_atom.GetIdx())
            amap = {}
            submol = Chem.PathToSubmol(mol, env, atomMap=amap)
            if submol.HasSubstructMatch(aminoethyl_pattern):
                attachment_found = True
                break
    if not attachment_found:
        return False, "Molecule does not have an aminoethyl side chain attached to the ring"

    return False, "Molecule does not match catecholamine criteria"