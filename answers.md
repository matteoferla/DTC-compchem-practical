See 2022 notebooks for the answers


```python
cen: npt.NDArray[np.float64] = get_centroid(mol.GetConformer())
shift_to_origin: npt.NDArray[np.float64] = create_translation_matrix(*(-cen))
shift_from_origin: npt.NDArray[np.float64] = create_translation_matrix(*(+cen))
rotation: npt.NDArray[np.float64] = create_rotation_matrix([1,0,0], angle=np.pi / 2)
# combine the affine transform matrices by matrix multiplication
_r: npt.NDArray[np.float64] = np.matmul(shift_to_origin, rotation)
rotation_on_spot: npt.NDArray[np.float64] = np.matmul(_r, shift_from_origin)

transform(mol.GetConformer(), rotation_on_spot )
```