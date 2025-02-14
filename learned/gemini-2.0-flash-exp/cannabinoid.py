"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    This function identifies characteristics like aromatic rings, long alkyl chains,
    oxygen atoms, and specific functional groups. It also tries to handle both
    classical and non-classical cannabinoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    score = 0
    reasons = []

    # 1. Check for dibenzopyran core (classical cannabinoids)
    dibenzopyran_pattern = Chem.MolFromSmarts("c1cc2c3c(cc1)Oc1ccccc1C3CC2")
    if mol.HasSubstructMatch(dibenzopyran_pattern):
        score += 3
        reasons.append("Contains dibenzopyran core")

    # 2. Check for other aromatic ring systems with oxygen (heterocycles)
    aromatic_oxygen_pattern1 = Chem.MolFromSmarts("c1ccc[o,n]c1")
    aromatic_oxygen_pattern2 = Chem.MolFromSmarts("c1cc[o,n]cc1")
    if mol.HasSubstructMatch(aromatic_oxygen_pattern1) or mol.HasSubstructMatch(aromatic_oxygen_pattern2):
         score += 2
         reasons.append("Contains aromatic heterocycle")

    # 3. Check for long alkyl chains (fatty acid or isoprenoid like) with at least one double bond
    # This pattern looks for a carbon chain of at least 4 carbons with at least one double bond, potentially with other heteroatoms.
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]=[CX3,CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    if len(fatty_acid_matches) > 0 :
        score += 2
        reasons.append("Contains fatty acid chain")

    # 4. Check for oxygen atoms
    has_oxygen_atoms = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if has_oxygen_atoms:
        score += 1
        reasons.append("Contains oxygen atoms")

    # 5. Check for key functional groups (esters, ethers, amides, alcohols, carbonyls, carboxyls, epoxides)

    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if mol.HasSubstructMatch(ester_pattern):
        score +=1
        reasons.append("Contains ester group")

    ether_pattern = Chem.MolFromSmarts("[OX2]-[CX4]")
    if mol.HasSubstructMatch(ether_pattern):
        score += 1
        reasons.append("Contains ether group")


    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX2]")
    if mol.HasSubstructMatch(amide_pattern):
        score += 1
        reasons.append("Contains amide group")

    alcohol_pattern = Chem.MolFromSmarts("[OX2][H]")
    if mol.HasSubstructMatch(alcohol_pattern):
        score +=1
        reasons.append("Contains alcohol group")

    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        score +=1
        reasons.append("Contains carbonyl group")
        
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O[H]")
    if mol.HasSubstructMatch(carboxyl_pattern):
        score += 1
        reasons.append("Contains carboxyl group")

    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(epoxide_pattern):
      score += 1
      reasons.append("Contains epoxide group")
    
    # Apply threshold
    if score >= 4:
       return True, ", ".join(reasons)
    else:
       return False, "Does not fit the criteria for a cannabinoid structure."