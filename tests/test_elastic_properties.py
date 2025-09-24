import pytest
from geomechpy.elastic_properties import ElasticProperties, ElasticPropertiesConverter

CONVERTER_TOLERANCE = 1e-2


def test_elastic_properties_converter_complex_approach() -> None:
    reference_model = ElasticProperties(
        bulk_modulus=40.56,
        youngs_modulus=32.55,
        lame_parameter=32.60,
        shear_modulus=11.91,
        poissons_ratio=0.366,
        p_wave_modulus=56.43,
    )

    calculated_models = [
        ElasticPropertiesConverter.from_bulk_and_youngs_modulus(reference_model.bulk_modulus, reference_model.youngs_modulus),
        ElasticPropertiesConverter.from_bulk_and_lame_modulus(reference_model.bulk_modulus, reference_model.lame_parameter),
        ElasticPropertiesConverter.from_bulk_and_shear_modulus(reference_model.bulk_modulus, reference_model.shear_modulus),
        ElasticPropertiesConverter.from_bulk_and_poissons_ratio_modulus(reference_model.bulk_modulus, reference_model.poissons_ratio),
        ElasticPropertiesConverter.from_bulk_and_p_wave_modulus(reference_model.bulk_modulus, reference_model.p_wave_modulus),
        ElasticPropertiesConverter.from_youngs_and_lame_modulus(reference_model.youngs_modulus, reference_model.lame_parameter),
        ElasticPropertiesConverter.from_youngs_and_shear_modulus(reference_model.youngs_modulus, reference_model.shear_modulus),
        ElasticPropertiesConverter.from_youngs_and_poissons_ratio_modulus(reference_model.youngs_modulus, reference_model.poissons_ratio),
        ElasticPropertiesConverter.from_youngs_and_p_wave_modulus(reference_model.youngs_modulus, reference_model.p_wave_modulus),
        ElasticPropertiesConverter.from_lame_and_shear_modulus(reference_model.lame_parameter, reference_model.shear_modulus),
        ElasticPropertiesConverter.from_lame_and_poissons_ratio_modulus(reference_model.lame_parameter, reference_model.poissons_ratio),
        ElasticPropertiesConverter.from_lame_and_p_wave_modulus(reference_model.lame_parameter, reference_model.p_wave_modulus),
        ElasticPropertiesConverter.from_shear_and_poissons_ratio_modulus(reference_model.shear_modulus, reference_model.poissons_ratio),
        ElasticPropertiesConverter.from_shear_and_p_wave_modulus(reference_model.shear_modulus, reference_model.p_wave_modulus),
        ElasticPropertiesConverter.from_poissons_ratio_and_p_wave_modulus(reference_model.poissons_ratio, reference_model.p_wave_modulus),
    ]

    for calculated_model in calculated_models:
        assert reference_model.bulk_modulus == pytest.approx(calculated_model.bulk_modulus, rel=CONVERTER_TOLERANCE)
        assert reference_model.youngs_modulus == pytest.approx(calculated_model.youngs_modulus, rel=CONVERTER_TOLERANCE)
        assert reference_model.lame_parameter == pytest.approx(calculated_model.lame_parameter, rel=CONVERTER_TOLERANCE)
        assert reference_model.shear_modulus == pytest.approx(calculated_model.shear_modulus, rel=CONVERTER_TOLERANCE)
        assert reference_model.poissons_ratio == pytest.approx(calculated_model.poissons_ratio, rel=CONVERTER_TOLERANCE)
        assert reference_model.p_wave_modulus == pytest.approx(calculated_model.p_wave_modulus, rel=CONVERTER_TOLERANCE)


def test_elastic_properties_convert_from_velosity() -> None:
    calculated_model = ElasticPropertiesConverter.from_velocity(30, 20, 2.5)

    assert 916.66 == pytest.approx(calculated_model.bulk_modulus, rel=CONVERTER_TOLERANCE)
    assert 2200.0 == pytest.approx(calculated_model.youngs_modulus, rel=CONVERTER_TOLERANCE)
    assert 250.0 == pytest.approx(calculated_model.lame_parameter, rel=CONVERTER_TOLERANCE)
    assert 1000.0 == pytest.approx(calculated_model.shear_modulus, rel=CONVERTER_TOLERANCE)
    assert 0.1 == pytest.approx(calculated_model.poissons_ratio, rel=CONVERTER_TOLERANCE)
    assert 2250.0 == pytest.approx(calculated_model.p_wave_modulus, rel=CONVERTER_TOLERANCE)


def test_elastic_properties_convert_from_slowness() -> None:
    calculated_model = ElasticPropertiesConverter.from_slowness(1000000, 2000000, 2500000)

    assert 154838.40 == pytest.approx(calculated_model.bulk_modulus, rel=CONVERTER_TOLERANCE)
    assert 154838.40 == pytest.approx(calculated_model.youngs_modulus, rel=CONVERTER_TOLERANCE)
    assert 116128.8 == pytest.approx(calculated_model.lame_parameter, rel=CONVERTER_TOLERANCE)
    assert 58064.40 == pytest.approx(calculated_model.shear_modulus, rel=CONVERTER_TOLERANCE)
    assert 0.33 == pytest.approx(calculated_model.poissons_ratio, rel=CONVERTER_TOLERANCE)
    assert 232257.6 == pytest.approx(calculated_model.p_wave_modulus, rel=CONVERTER_TOLERANCE)
