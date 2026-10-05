from urllib.parse import quote

from pancgi_contract import integer


def representative_allele_id(member):
    names = []
    for field in ('hal_genome', 'hal_sequence'):
        value = member[field]
        if not isinstance(value, str) or not value or value != value.strip() or '\x00' in value:
            raise ValueError(f'Invalid representative {field}')
        names.append(quote(value, safe='-._~', encoding='utf-8', errors='strict'))
    start = integer(str(member['start0']), 'Representative start0')
    end = integer(str(member['end0']), 'Representative end0')
    if end <= start:
        raise ValueError('Representative CGI interval must have positive length')
    return f'{names[0]}_{names[1]}:{start}:{end}'


def allele_id_map(catalogue, members):
    result, used = {}, set()
    for row in catalogue:
        internal = row['allele_id']
        if not internal or internal in result:
            raise ValueError('Missing or duplicate internal allele identifier')
        representative = row['allele_rep_fid']
        if representative not in members:
            raise ValueError('Unknown allele representative member')
        public = representative_allele_id(members[representative])
        if public in used:
            raise ValueError(f'Public allele identifier collision: {public}')
        result[internal] = public
        used.add(public)
    return result
