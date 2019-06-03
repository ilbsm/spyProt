#!/usr/bin/env python
import urllib.request, urllib.error, urllib.parse


def fetchReleasedPdbs(from_date, to_date=''):
    '''
        Download list of pdb ids released in range from_date - to_date (or in day from_date). Return False if empty list
    '''
    if to_date == '':
        to_date = from_date
    '''
    date in format 2009-07-01
    '''
    url = 'http://www.rcsb.org/pdb/rest/search'
    queryText = """
<orgPdbCompositeQuery version="1.0">
<queryRefinement>
<queryRefinementLevel>0</queryRefinementLevel>
<orgPdbQuery>

<queryType>org.pdb.query.simple.ReleaseDateQuery</queryType>
<pdbx_audit_revision_history.revision_date.comparator>between</pdbx_audit_revision_history.revision_date.comparator>
<pdbx_audit_revision_history.revision_date.min>%s</pdbx_audit_revision_history.revision_date.min>
<pdbx_audit_revision_history.revision_date.max>%s</pdbx_audit_revision_history.revision_date.max>
<pdbx_audit_revision_history.ordinal.comparator>=</pdbx_audit_revision_history.ordinal.comparator>
<pdbx_audit_revision_history.ordinal.value>1</pdbx_audit_revision_history.ordinal.value>
</orgPdbQuery>
</queryRefinement>
<queryRefinement>
<queryRefinementLevel>1</queryRefinementLevel>
<conjunctionType>and</conjunctionType>
<orgPdbQuery>
<queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
<description>Chain Type: there is a Protein chain</description>
<containsProtein>Y</containsProtein>
<containsDna>?</containsDna>
<containsRna>?</containsRna>
<containsHybrid>?</containsHybrid>
</orgPdbQuery>
</queryRefinement>
</orgPdbCompositeQuery>
    """ % (from_date, to_date)

    req = urllib.request.Request(url, data=queryText.encode("utf-8"), method='POST')
    f = urllib.request.urlopen(req)
    result = [i.strip().decode("utf-8") for i in f.readlines()]
    f.close()
    if result:
        return result
    else:
        return False


if __name__ == "__main__":
    print(fetchReleasedPdbs("2017-10-01"))
#    print fetchReleasedPdbs("2017-08-1","2017-08-20") 

