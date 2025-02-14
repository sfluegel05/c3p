"""
Classifies: CHEBI:33839 macromolecule
"""
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is defined as a molecule of high relative molecular mass,
    with its structure comprising the multiple repetition of units derived from
    smaller molecules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Check for common macromolecular repeat units
    repeat_unit_patterns = ['C(NC([C@H](NC([C@H](NC([C@H](NC([C@H](NC(C(NC([C@H](NC([C@H](NC([C@H](NC([C@H](NC(C)=O)C)=O)C)CC)C(=O)N)=O)CC)=O)CCC)=O)CC)=O)CC)=O)CC)=O)CC',  # peptide
                            '[OC@H]([C@H]([C@H]([C@H](O)O)O)O)[C@@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O',  # polysaccharide
                            'C(C(C(=O)OC)C)C(=O)OC',  # polyester
                            'C(CC(C)C)C(=O)N[C@@H](CSSC[C@H](C(=O)O)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]([C@@H](C)CC)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](CO)NC(=O)[C@@H](CCC(=O)O)NC(=O)[C@@H](CC(=O)O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@@H](CO)NC(=O)[C@@H](C(C)C)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@@H]([C@@H](C)CC)NC(=O)[C@H](CC(N)=O)NC(=O)[C@@H](CCCCN)NC(=O)[C@@H](CO)NC(=O)[C@@H](CO)NC(=O)QC(=O)NC(=O)CCC(=O)O)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N1CCCC[C@H](NC(=O)[C@@H](N)CC(C)C)C(=O)NCC(=O)N']  # non-ribosomal peptide

    for pattern in repeat_unit_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Contains repeat unit pattern: {pattern}"

    # Check for long carbon chains
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if len(carbon_chain_matches) > 0:
        return True, "Contains long carbon chains"

    # Check molecular weight - macromolecules typically > 1000 Da
    if mol_wt > 1000:
        return True, f"Molecular weight of {mol_wt:.2f} Da suggests macromolecule"

    # Check for other common macromolecular features
    if re.search(r'(\([^\)]+\)\1\1)+', smiles):  # Repeating units in SMILES
        return True, "SMILES string contains repeating units"

    if len(smiles) > 200:  # Long SMILES strings often indicate macromolecules
        return True, "SMILES string is very long, indicating a macromolecule"

    return False, "No macromolecular features detected"