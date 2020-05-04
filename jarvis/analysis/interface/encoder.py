import os
import json
import datetime

from hashlib import sha1


class MontyEncoder(json.JSONEncoder):
    """
    A Json Encoder which supports the MSONable API, plus adds support for
    numpy arrays, datetime objects, bson ObjectIds (requires bson).
    Usage::
        # Add it as a *cls* keyword when using json.dump
        json.dumps(object, cls=MontyEncoder)
    """

    def default(self, o) -> dict:  # pylint: disable=E0202
        """
        Overriding default method for JSON encoding. This method does two
        things: (a) If an object has a to_dict property, return the to_dict
        output. (b) If the @module and @class keys are not in the to_dict,
        add them to the output automatically. If the object has no to_dict
        property, the default Python json encoder default method is called.
        Args:
            o: Python object.
        Return:
            Python dict representation.
        """
        if isinstance(o, datetime.datetime):
            return {"@module": "datetime", "@class": "datetime", "string": o.__str__()}
        if np is not None:
            if isinstance(o, np.ndarray):
                return {
                    "@module": "numpy",
                    "@class": "array",
                    "dtype": o.dtype.__str__(),
                    "data": o.tolist(),
                }
            if isinstance(o, np.generic):
                return o.item()
        if bson is not None:
            if isinstance(o, bson.objectid.ObjectId):
                return {"@module": "bson.objectid", "@class": "ObjectId", "oid": str(o)}

        try:
            d = o.as_dict()
            if "@module" not in d:
                d["@module"] = u"{}".format(o.__class__.__module__)
            if "@class" not in d:
                d["@class"] = u"{}".format(o.__class__.__name__)
            if "@version" not in d:
                try:
                    parent_module = o.__class__.__module__.split(".")[0]
                    module_version = import_module(
                        parent_module
                    ).__version__  # type: ignore
                    d["@version"] = u"{}".format(module_version)
                except (AttributeError, ImportError):
                    d["@version"] = None
            return d
        except AttributeError:
            return json.JSONEncoder.default(self, o)


class MontyDecoder(json.JSONDecoder):
    """
    A Json Decoder which supports the MSONable API. By default, the
    decoder attempts to find a module and name associated with a dict. If
    found, the decoder will generate a Pymatgen as a priority.  If that fails,
    the original decoded dictionary from the string is returned. Note that
    nested lists and dicts containing pymatgen object will be decoded correctly
    as well.
    Usage:
        # Add it as a *cls* keyword when using json.load
        json.loads(json_string, cls=MontyDecoder)
    """

    def process_decoded(self, d):
        """
        Recursive method to support decoding dicts and lists containing
        pymatgen objects.
        """
        if isinstance(d, dict):
            if "@module" in d and "@class" in d:
                modname = d["@module"]
                classname = d["@class"]
                if classname in MSONable.REDIRECT.get(modname, {}):
                    modname = MSONable.REDIRECT[modname][classname]["@module"]
                    classname = MSONable.REDIRECT[modname][classname]["@class"]
            else:
                modname = None
                classname = None
            if modname and modname not in ["bson.objectid", "numpy"]:
                if modname == "datetime" and classname == "datetime":
                    try:
                        dt = datetime.datetime.strptime(
                            d["string"], "%Y-%m-%d %H:%M:%S.%f"
                        )
                    except ValueError:
                        dt = datetime.datetime.strptime(
                            d["string"], "%Y-%m-%d %H:%M:%S"
                        )
                    return dt

                mod = __import__(modname, globals(), locals(), [classname], 0)
                if hasattr(mod, classname):
                    cls_ = getattr(mod, classname)
                    data = {k: v for k, v in d.items() if not k.startswith("@")}
                    if hasattr(cls_, "from_dict"):
                        return cls_.from_dict(data)
            elif np is not None and modname == "numpy" and classname == "array":
                return np.array(d["data"], dtype=d["dtype"])

            elif (
                (bson is not None)
                and modname == "bson.objectid"
                and classname == "ObjectId"
            ):
                return bson.objectid.ObjectId(d["oid"])

            return {
                self.process_decoded(k): self.process_decoded(v) for k, v in d.items()
            }

        if isinstance(d, list):
            return [self.process_decoded(x) for x in d]

        return d

    def decode(self, s):
        """
        Overrides decode from JSONDecoder.
        :param s: string
        :return: Object.
        """
        d = json.JSONDecoder.decode(self, s)
        return self.process_decoded(d)
