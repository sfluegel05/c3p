"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
But this code might not handle the complexity of glycosaminoglycans properly. Let's think about how to better identify polysaccharide structures and aminomonosaccharides.

Possible issues:
1. The code checks for multiple rings but doesn't confirm they are connected via glycosidic bonds.
2. The amino group detection may not be specific enough (could pick up non-saccharide amines).
3. The 20% threshold is arbitrary and may not align with the definition.

Improvements:
- Use more specific patterns for saccharide units (e.g., six-membered rings with oxygen).
- Check for glycosidic linkages (O connecting two rings).
- Ensure amino groups are part of the saccharide rings.
- Consider molecular weight to confirm it's a polysaccharide.

Let's try a better approach: