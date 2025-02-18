"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:64211 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by oxygen-containing groups in specific arrangements,
    including endocannabinoids (ethanolamides, glycerol esters) and phytocannabinoids (dibenzopyran derivatives).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for ethanolamide group (N-linked to carbonyl and ethanol)
    ethanolamide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H][CX4H2][CX4H2]O")
    if mol.HasSubstructMatch(ethanolamide_pattern):
        return True, "Contains ethanolamide group characteristic of endocannabinoids"

    # Check for glycerol-linked structures (ester/ether with long chain)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4](O[!H0])[CH2X4]")
    if mol.HasSubstructMatch(glycerol_pattern):
        # Verify at least one long chain (>=16 heavy atoms in connected chain)
        for match in mol.GetSubstructMatches(glycerol_pattern):
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, 4, match[1])
            atoms = {a for a in env if mol.GetAtomWithIdx(a).GetAtomicNum() == 6}
            if len(atoms) >= 14:
                return True, "Contains glycerol backbone with long-chain modification"

    # Check for dibenzopyran core (THC-like)
    dibenzopyran_core = Chem.MolFromSmarts("[C]1([C]2=[C]([C@@H]3[C@@H](CC[C@H]3O)C)[CH]=[CH][C]=2)[CH2][CH2][CH2]1")
    if mol.HasSubstructMatch(dibenzopyran_core):
        return True, "Contains dibenzopyran cannabinoid core"

    # Check for synthetic cannabinoid scaffolds (indole/indazole + carbonyl linker)
    synth_pattern = Chem.MolFromSmarts("[c]1[c][c][c](C(=O)[!C;!N])[c]([!O])[nH]1")  # Indole carbonyl
    if mol.HasSubstructMatch(synth_pattern):
        return True, "Contains synthetic cannabinoid scaffold"

    # Oxygen in heterocyclic ring (corrected from atomic number 16 to 8)
    oxygen_ring = any(atom.GetAtomicNum() == 8 and atom.IsInRing() for atom in mol.GetAtoms())
    if oxygen_ring:
        # Require minimum ring size of 5 and adjacent carbons
        rings = mol.GetRingInfo().AtomRings()
        for ring in rings:
            if any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring) and len(ring) >= 5:
                return True, "Contains oxygen in heterocyclic ring system"

    # Check for long unsaturated chain (>=18 carbons) with oxygen groups
    long_chain = Chem.MolFromSmarts("[CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4]")
    if mol.HasSubstructMatch(long_chain):
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        if oxygen_count >= 2:
            return True, "Long unsaturated chain with oxygen functional groups"

    return False, "Does not match known cannabinoid structural patterns"