import pathlib
from unittest.mock import patch

import numpy as np
import pytest

from lattice_mc import init_lattice

# TODO write integration tests for creating lattices


class TestInitLatticeDimensionality:
    """Tests for the dimensionality of lattices built by the initialisation routines"""

    def test_square_lattice_is_two_dimensional(self):
        assert init_lattice.square_lattice(3, 4, 1.0).dimensionality == 2

    def test_honeycomb_lattice_is_two_dimensional(self):
        assert init_lattice.honeycomb_lattice(2, 2, 1.0).dimensionality == 2

    def test_cubic_lattice_is_three_dimensional(self):
        assert init_lattice.cubic_lattice(3, 3, 3, 1.0).dimensionality == 3


class TestInitLatticeSizes:
    """Tests for the smallest numbers of repeat units the builders accept"""

    @pytest.mark.parametrize(
        ("build", "argument"),
        [
            (lambda: init_lattice.square_lattice(2, 3, 1.0), "a"),
            (lambda: init_lattice.square_lattice(3, 1, 1.0), "b"),
            (lambda: init_lattice.cubic_lattice(2, 3, 3, 1.0), "a"),
            (lambda: init_lattice.cubic_lattice(3, 2, 3, 1.0), "b"),
            (lambda: init_lattice.cubic_lattice(3, 3, 1, 1.0), "c"),
            (lambda: init_lattice.honeycomb_lattice(1, 2, 1.0), "a"),
            (lambda: init_lattice.honeycomb_lattice(2, 0, 1.0), "b"),
        ],
        ids=["square a", "square b", "cubic a", "cubic b", "cubic c", "honeycomb a", "honeycomb b"],
    )
    def test_rejects_too_few_repeat_units(self, build, argument):
        with pytest.raises(ValueError, match=f"^{argument} must be at least"):
            build()

    @pytest.mark.parametrize(
        ("build", "number_of_sites"),
        [
            (lambda: init_lattice.square_lattice(3, 3, 1.0), 9),
            (lambda: init_lattice.cubic_lattice(3, 3, 3, 1.0), 27),
            (lambda: init_lattice.honeycomb_lattice(2, 1, 1.0), 8),
        ],
        ids=["square", "cubic", "honeycomb"],
    )
    def test_builds_the_smallest_accepted_lattice(self, build, number_of_sites):
        assert build().number_of_sites == number_of_sites


class TestInitLattice:
    """Test for Specific Lattice initialisation routines"""

    @patch("lattice_mc.lattice_site.Site")
    @patch("lattice_mc.lattice.Lattice")
    def test_square_lattice(self, mock_lattice, mock_site):
        mock_site.side_effect = range(3 * 4)
        mock_lattice.return_value = "foo"
        a = 3
        b = 4
        spacing = 1.0
        lattice = init_lattice.square_lattice(a, b, spacing)
        expected_site_calls = [
            [1, np.array([0.0, 0.0, 0.0]), [3, 2, 10, 4], 0.0, "L"],
            [2, np.array([1.0, 0.0, 0.0]), [1, 3, 11, 5], 0.0, "L"],
            [3, np.array([2.0, 0.0, 0.0]), [2, 1, 12, 6], 0.0, "L"],
            [4, np.array([0.0, 1.0, 0.0]), [6, 5, 1, 7], 0.0, "L"],
            [5, np.array([1.0, 1.0, 0.0]), [4, 6, 2, 8], 0.0, "L"],
            [6, np.array([2.0, 1.0, 0.0]), [5, 4, 3, 9], 0.0, "L"],
            [7, np.array([0.0, 2.0, 0.0]), [9, 8, 4, 10], 0.0, "L"],
            [8, np.array([1.0, 2.0, 0.0]), [7, 9, 5, 11], 0.0, "L"],
            [9, np.array([2.0, 2.0, 0.0]), [8, 7, 6, 12], 0.0, "L"],
            [10, np.array([0.0, 3.0, 0.0]), [12, 11, 7, 1], 0.0, "L"],
            [11, np.array([1.0, 3.0, 0.0]), [10, 12, 8, 2], 0.0, "L"],
            [12, np.array([2.0, 3.0, 0.0]), [11, 10, 9, 3], 0.0, "L"],
        ]
        for c, e in zip(mock_site.mock_calls, expected_site_calls, strict=True):
            assert c[1][0] == e[0]  # site number
            np.testing.assert_array_equal(c[1][1], e[1])  # site coordinates
            assert c[1][2] == e[2]  # neighbour lists
            assert c[1][3] == e[3]  # ?
            assert c[1][4] == e[4]  # site label
        assert mock_lattice.mock_calls[0][1][0] == list(range(3 * 4))
        assert lattice == "foo"
        np.testing.assert_array_equal(mock_lattice.mock_calls[0][2]["cell_lengths"], np.array([3.0, 4.0, 0.0]))

    @patch("lattice_mc.lattice_site.Site")
    @patch("lattice_mc.lattice.Lattice")
    def test_honeycomb_lattice(self, mock_lattice, mock_site):
        a, b = 2, 1
        mock_site.side_effect = range(a * b * 4)
        mock_lattice.return_value = "foo"
        spacing = 1.0
        lattice = init_lattice.honeycomb_lattice(a, b, spacing)
        expected_site_calls = [
            [1, np.array([0.0, 0.0, 0.0]), [2, 4, 6], 0.0, "L"],
            [2, np.array([0.8660254, 0.5, 0.0]), [1, 3, 5], 0.0, "L"],
            [3, np.array([0.8660254, 1.5, 0.0]), [2, 4, 8], 0.0, "L"],
            [4, np.array([0.0, 2.0, 0.0]), [1, 3, 7], 0.0, "L"],
            [5, np.array([1.73205081, 0.0, 0.0]), [2, 6, 8], 0.0, "L"],
            [6, np.array([2.59807621, 0.5, 0.0]), [5, 1, 7], 0.0, "L"],
            [7, np.array([2.59807621, 1.5, 0.0]), [8, 6, 4], 0.0, "L"],
            [8, np.array([1.73205081, 2.0, 0.0]), [3, 5, 7], 0.0, "L"],
        ]
        for c, e in zip(mock_site.mock_calls, expected_site_calls, strict=True):
            assert c[1][0] == e[0]  # site number
            np.testing.assert_array_almost_equal(c[1][1], e[1])  # site coordinates
            assert sorted(c[1][2]) == sorted(e[2])  # neighbour lists
            assert c[1][3] == e[3]  # ?
            assert c[1][4] == e[4]  # site label
        assert mock_lattice.mock_calls[0][1][0] == list(range(2 * 1 * 4))
        assert lattice == "foo"
        np.testing.assert_array_almost_equal(
            mock_lattice.mock_calls[0][2]["cell_lengths"], np.array([3.464102, 3.0, 0.0])
        )

    @patch("lattice_mc.lattice_site.Site")
    @patch("lattice_mc.lattice.Lattice")
    def test_honeycomb_lattice_periodic_boundaries(self, mock_lattice, mock_site):
        a, b = 3, 3
        mock_site.side_effect = range(a * b * 4)
        mock_lattice.return_value = "foo"
        spacing = 1.0
        init_lattice.honeycomb_lattice(a, b, spacing)
        expected_neighbours = {1: [2, 26, 12], 12: [1, 35, 11], 35: [36, 34, 12], 26: [27, 25, 1]}
        for c in mock_site.mock_calls:
            site_number = c[1][0]
            if site_number in expected_neighbours:
                site_neighbours = c[1][2]
                assert sorted(site_neighbours) == sorted(expected_neighbours[site_number])

    @patch("lattice_mc.lattice_site.Site")
    @patch("lattice_mc.lattice.Lattice")
    def test_honeycomb_lattice_alternating_sites(self, mock_lattice, mock_site):
        a, b = 3, 3
        mock_site.side_effect = range(a * b * 4)
        mock_lattice.return_value = "foo"
        spacing = 1.0
        init_lattice.honeycomb_lattice(a, b, spacing, alternating_sites=True)
        for c in mock_site.mock_calls:
            site_label = c[1][4]
            site_neighbours = c[1][2]
            for n in site_neighbours:
                neighbour_label = mock_site.mock_calls[n - 1][1][4]
                assert set([site_label, neighbour_label]) == set(["A", "B"])

    @patch("lattice_mc.lattice_site.Site")
    @patch("lattice_mc.lattice.Lattice")
    def test_cubic_lattice(self, mock_lattice, mock_site):
        a, b, c = 3, 3, 3
        spacing = 1.0
        mock_site.side_effect = range(3 * 3 * 3)
        mock_lattice.return_value = "bar"
        lattice = init_lattice.cubic_lattice(a, b, c, spacing)

        def number(x, y, z):
            # Sites are numbered from 1, with x varying fastest.
            return 1 + x % 3 + 3 * (y % 3) + 9 * (z % 3)

        expected_site_calls = [
            [
                number(x, y, z),
                np.array([x, y, z], dtype=float),
                [
                    number(x - 1, y, z),
                    number(x + 1, y, z),
                    number(x, y - 1, z),
                    number(x, y + 1, z),
                    number(x, y, z - 1),
                    number(x, y, z + 1),
                ],
                0.0,
                "L",
            ]
            for z in range(3)
            for y in range(3)
            for x in range(3)
        ]
        for call, e in zip(mock_site.mock_calls, expected_site_calls, strict=True):
            assert call[1][0] == e[0]  # site number
            np.testing.assert_array_equal(call[1][1], e[1])  # site coordinates
            assert call[1][2] == e[2]  # neighbour lists
            assert call[1][3] == e[3]  # ?
            assert call[1][4] == e[4]  # site label
        assert mock_lattice.mock_calls[0][1][0] == list(range(3 * 3 * 3))
        assert lattice == "bar"
        np.testing.assert_array_equal(mock_lattice.mock_calls[0][2]["cell_lengths"], np.array([a, b, c]))


RING = """3

site: 1
centre: 0.0 0.0 0.0
neighbours: 2 3
label: A

site: 2
centre: 1.0 0.0 0.0
neighbours: 1 3
label: A

site: 3
centre: 2.0 0.0 0.0
neighbours: 1 2
label: A
"""


def write_sites_file(tmp_path, text):
    path = tmp_path / "sites.dat"
    path.write_text(text)
    return str(path)


def load(tmp_path, text):
    return init_lattice.lattice_from_sites_file(write_sites_file(tmp_path, text), cell_lengths=[3.0, 10.0, 10.0])


class TestLatticeFromSitesFile:
    def test_reads_sites(self, tmp_path):
        lattice = load(tmp_path, RING)
        site = lattice.site_with_id(2)
        np.testing.assert_array_equal(site.r, [1.0, 0.0, 0.0])
        assert site.neighbours == [1, 3]
        assert site.label == "A"
        assert site.energy == 0.0
        np.testing.assert_array_equal(lattice.cell_lengths, [3.0, 10.0, 10.0])

    def test_reads_energy(self, tmp_path):
        lattice = load(tmp_path, RING.replace("label: A\n", "label: A\nenergy: -0.5\n"))
        assert {site.energy for site in lattice.sites} == {-0.5}

    def test_reads_energy_in_exponent_notation(self, tmp_path):
        lattice = load(tmp_path, RING.replace("label: A\n", "label: A\nenergy: 1e-3\n"))
        assert {site.energy for site in lattice.sites} == {1e-3}

    def test_allows_extra_whitespace_after_a_key(self, tmp_path):
        lattice = load(tmp_path, RING.replace("label: A\n", "label:   A\nenergy:  0.1\n"))
        assert {(site.label, site.energy) for site in lattice.sites} == {("A", 0.1)}

    def test_allows_blank_lines_containing_whitespace(self, tmp_path):
        lattice = load(tmp_path, RING.replace("\n\n", "\n   \n"))
        assert lattice.number_of_sites == 3

    def test_accepts_american_spellings_and_any_line_order(self, tmp_path):
        text = RING.replace("centre:", "center:").replace("neighbours:", "neighbors:")
        text = text.replace(
            "site: 1\ncenter: 0.0 0.0 0.0\nneighbors: 2 3\nlabel: A",
            "label: A\nneighbors: 2 3\ncenter: 0.0 0.0 0.0\nsite: 1",
        )
        lattice = load(tmp_path, text)
        assert lattice.site_with_id(1).neighbours == [2, 3]

    def test_ignores_legacy_vertices_lines(self, tmp_path):
        lattice = load(tmp_path, RING.replace("label: A\n", "label: A\nvertices: 1 315 435 649\n"))
        assert lattice.number_of_sites == 3

    @pytest.mark.parametrize(
        ("old", "new", "message"),
        [
            ("3\n", "4\n", "4 sites"),
            ("3\n", "three\n", "number of sites"),
            ("site: 3", "site: 2", "site 2"),
            ("label: A\n\nsite: 2", "colour: red\nlabel: A\n\nsite: 2", "colour"),
            ("label: A\n\nsite: 2", "label: A\nlabel: A\n\nsite: 2", "'label' appears more than once"),
            ("label: A\n\nsite: 2", "label: A\nnonsense\n\nsite: 2", "key: value"),
            ("neighbours: 1 3\nlabel: A\n", "neighbours: 1 3\n", "label"),
            ("centre: 1.0 0.0 0.0", "centre: 1.0 0.0", "three coordinates"),
            ("centre: 2.0 0.0 0.0", "centre: 2.0 x 0.0", "Could not read"),
            ("label: A\n", "label: A B\n", "single word"),
            ("neighbours: 1 2\n", "neighbours: 1 2 4\n", "site 4"),
            ("neighbours: 2 3\n", "neighbours: 2\n", "site 3"),
            ("neighbours: 2 3\n", "neighbours: 1 2 3\n", "itself"),
            ("neighbours: 2 3\n", "neighbours: 2 2 3\n", "site 2 more than once"),
            ("neighbours: 1 3\nlabel: A\n", "neighbours: 1 3\nlabel: A\nenergy: nan\n", "site 2 must be finite"),
        ],
        ids=[
            "count does not match",
            "header is not an integer",
            "duplicate site number",
            "unknown key",
            "repeated key",
            "line without a colon",
            "missing label",
            "centre without three coordinates",
            "non-numeric value",
            "multi-word label",
            "neighbour does not exist",
            "one-way neighbours",
            "self-neighbour",
            "repeated neighbour",
            "non-finite energy",
        ],
    )
    def test_raises_ValueError_for_invalid_file(self, tmp_path, old, new, message):
        assert old in RING
        with pytest.raises(ValueError, match=message):
            load(tmp_path, RING.replace(old, new, 1))

    def test_reads_different_energies_for_sites_with_one_label(self, tmp_path):
        text = RING.replace(
            "site: 1\ncentre: 0.0 0.0 0.0\nneighbours: 2 3\nlabel: A\n",
            "site: 1\ncentre: 0.0 0.0 0.0\nneighbours: 2 3\nlabel: A\nenergy: 0.2\n",
        )
        lattice = load(tmp_path, text)
        assert [lattice.site_with_id(n).energy for n in (1, 2, 3)] == [0.2, 0.0, 0.0]

    def test_reads_the_example_llzo_sites_file(self):
        path = str(pathlib.Path(__file__).resolve().parents[1] / "examples" / "llzo_lattice_site_list.dat")
        lattice = init_lattice.lattice_from_sites_file(path, cell_lengths=[49.0672361] * 3)
        assert lattice.number_of_sites == 576
        assert lattice.site_labels == {"O", "T"}
