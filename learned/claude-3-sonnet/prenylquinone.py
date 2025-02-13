"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: CHEBI:51829 prenylquinone

A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.
They are characterized by a quinone core (benzoquinone, naphthoquinone, anthraquinone, etc.)
with one or more prenyl/polyprenyl chains attached.
"""

from rdkit import Chem
from rdkit.Chem import rdqueries, rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for quinone core
    quinone_patterns = ['*1=,=,*2*3=,=,*4*=1*=2*=3*=4', # Benzoquinone
                        'c1c(=O)c2ccccc2c(=O)c1', # Naphthoquinone
                        'c1c(=O)c2ccccc2c(=O)c3ccccc13'] # Anthraquinone
    quinone_match = any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in quinone_patterns)
    if not quinone_match:
        return False, "No quinone core found"

    # Check for prenyl/polyprenyl chains
    prenyl_pattern = rdqueries.PolyprenylQueries().HavingPolyprenylSidechains()
    prenyl_match = mol.HasSubstructMatch(prenyl_pattern)
    if not prenyl_match:
        return False, "No prenyl/polyprenyl chains found"

    # Check if prenyl chains are attached to quinone core
    quinone_atoms = set(hit for p in quinone_patterns
                            for hit in mol.GetSubstructMatches(Chem.MolFromSmarts(p)))
    prenyl_atoms = set(hit for hit in mol.GetSubstructMatches(prenyl_pattern))
    if not prenyl_atoms.intersection(quinone_atoms):
        return False, "Prenyl chains not attached to quinone core"

    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1500:
        return False, "Molecular weight outside expected range for prenylquinones"

    return True, "Contains a quinone core with prenyl/polyprenyl substituents"