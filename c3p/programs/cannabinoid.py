"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by specific structural features such as
    a resorcinol moiety with an alkyl side chain, dibenzopyran ring system,
    or long-chain polyunsaturated fatty acid derivatives linked to ethanolamine or glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for phytocannabinoids (e.g., THC, CBD)
    # Resorcinol moiety with alkyl side chain
    resorcinol_alkyl_pattern = Chem.MolFromSmarts('Oc1cc(O)ccc1CCCC')
    # Dibenzopyran ring system
    dibenzopyran_pattern = Chem.MolFromSmarts('c1cc2ccc(O)cc2oc1')

    # Define patterns for endocannabinoids (e.g., anandamide, 2-AG)
    # Long-chain fatty acid amide linked to ethanolamine
    fatty_acid_amide_pattern = Chem.MolFromSmarts('C(=O)NCCO')
    # Long-chain fatty acid ester linked to glycerol
    fatty_acid_glycerol_pattern = Chem.MolFromSmarts('C(=O)OCC(O)CO')

    # Check for phytocannabinoid features
    if mol.HasSubstructMatch(resorcinol_alkyl_pattern):
        return True, "Contains resorcinol moiety with alkyl side chain characteristic of phytocannabinoids"

    if mol.HasSubstructMatch(dibenzopyran_pattern):
        return True, "Contains dibenzopyran ring system characteristic of phytocannabinoids"

    # Check for endocannabinoid features
    if mol.HasSubstructMatch(fatty_acid_amide_pattern):
        # Check for long-chain fatty acid (at least 18 carbons)
        chain_lengths = []
        for amide_match in mol.GetSubstructMatches(fatty_acid_amide_pattern):
            carbon_chain = []
            amide_carbon = amide_match[0]
            atom = mol.GetAtomWithIdx(amide_carbon)
            # Traverse the carbon chain
            while True:
                neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() not in carbon_chain]
                if neighbors:
                    atom = neighbors[0]
                    carbon_chain.append(atom.GetIdx())
                else:
                    break
            chain_lengths.append(len(carbon_chain))
        if chain_lengths and max(chain_lengths) >= 16:
            return True, "Contains long-chain fatty acid amide characteristic of endocannabinoids"
        else:
            return False, "Amide group found but fatty acid chain too short"

    if mol.HasSubstructMatch(fatty_acid_glycerol_pattern):
        # Check for long-chain fatty acid (at least 18 carbons)
        chain_lengths = []
        for ester_match in mol.GetSubstructMatches(fatty_acid_glycerol_pattern):
            carbon_chain = []
            ester_carbon = ester_match[0]
            atom = mol.GetAtomWithIdx(ester_carbon)
            # Traverse the carbon chain
            while True:
                neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() not in carbon_chain]
                if neighbors:
                    atom = neighbors[0]
                    carbon_chain.append(atom.GetIdx())
                else:
                    break
            chain_lengths.append(len(carbon_chain))
        if chain_lengths and max(chain_lengths) >= 16:
            return True, "Contains long-chain fatty acid ester characteristic of endocannabinoids"
        else:
            return False, "Ester group found but fatty acid chain too short"

    return False, "No characteristic cannabinoid structural features found"