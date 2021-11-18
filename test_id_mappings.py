from spyprot import generate_mappings, PDB_Uniprot


def test_mapping():
    generate_mappings(just_pdb_to_uniprot=False)
    res = PDB_Uniprot("1j85")
    assert res == [('1j85 A', 'P44868')]
    result = PDB_Uniprot("P10970")
    assert str(result) == "['1ffk O', '1jj2 Q', '1k73 S', '1k8a S', '1k9m S', '1kc8 S', '1kd1 S', '1kqs Q', '1m1k S', '1m90 S', '1n8r S', '1nji S', '1q7y S', '1q81 S', '1q82 S', '1q86 S', '1qvf Q', '1qvg Q', '1s72 R', '1vq4 R', '1vq5 R', '1vq6 R', '1vq7 R', '1vq8 R', '1vq9 R', '1vqk R', '1vql R', '1vqm R', '1vqn R', '1vqo R', '1vqp R', '1w2b Q', '1yhq R', '1yi2 R', '1yij R', '1yit R', '1yj9 R', '1yjn R', '1yjw R', '2otj R', '2otl R', '2qa4 R', '2qex R', '3cc2 R', '3cc4 R', '3cc7 R', '3cce R', '3ccj R', '3ccl R', '3ccm R', '3ccq R', '3ccr R', '3ccs R', '3ccu R', '3ccv R', '3cd6 R', '3cma R', '3cme R', '3cpw Q', '3cxc Q', '3g4s R', '3g6e R', '3g71 R', '3i55 R', '3i56 R', '3ow2 Q', '4adx R', '4v9f R']"
