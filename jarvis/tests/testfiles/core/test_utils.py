from jarvis.core.utils import recast_array_on_uniq_array_elements,stringdict_to_xml,array_to_string
def test_utils():
   info={'a':'b','c':'d'}
   sd=stringdict_to_xml(info)
   sd=stringdict_to_xml(info,enforce_string=True)
   ar=[1,2,3,4,5]
   sarr=array_to_string(ar)
