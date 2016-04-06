from xml.dom.minidom import parse
import json

"""
A quick and dirty Python 3 script to map SBML reaction and compound IDs
to their KEGG counterparts.
Prints the map in JSON format.
"""


def sbml_parse(filename):
    """

    :param filename: the path to a valid SBML file.
    :return: a dictionary of the form model_id -> KEGG_id
    """

    file = open(filename)
    xml_dom = parse(file)

    # Find XML nodes
    species = xml_dom.getElementsByTagName("species")
    reactions = xml_dom.getElementsByTagName("reaction")
    node_list = species + reactions

    d = {}  # Empty dictionary of IDs
    for node in node_list:
        model_id = node.getAttribute("id")
        resources = node.getElementsByTagName("rdf:li")
        for r in resources:
            path = r.getAttribute("rdf:resource")
            # Filter out things that are not compounds or reactions:
            if "kegg.compound" in path or "kegg.reaction" in path:
                kegg_id = path.split("/")[-1]
                d[model_id] = kegg_id
    return d


if __name__ == "__main__":
    d = sbml_parse("yeast_7.00.xml")
    # To get a KEGG_id -> model_id map, inverse the dict:
    # d = {v: k for k, v in d.items()}
    print(json.dumps(d, indent=4, sort_keys=True))