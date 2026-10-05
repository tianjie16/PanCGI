MECHANISMS = ('INS-internal', 'INS-junction', 'INS-spanning', 'DEL-junction', 'INS+DEL', 'non-SV')
RELATIONSHIPS = {'inside_ins': 'INS-internal', 'partial_ins': 'INS-junction',
                 'contains_ins': 'INS-spanning', 'contains_del': 'DEL-junction'}


def mechanism(events, role='sample'):
    if role != 'sample':
        if events:
            raise ValueError('Reference-role member has unexpected SV events')
        return dict(mechanism_group='', mechanism_class='', mechanism_status='not_applicable_reference_role')
    classes = set()
    for event in events:
        relation = event['class']
        if relation not in RELATIONSHIPS:
            raise ValueError(f'Unsupported event relationship: {relation}')
        expected_type = 'DEL' if relation == 'contains_del' else 'INS'
        if event['svtype'] != expected_type:
            raise ValueError('Event type and relationship disagree')
        classes.add(RELATIONSHIPS[relation])
    ins = next((c for c in MECHANISMS[:3] if c in classes), '')
    result = 'INS+DEL' if ins and 'DEL-junction' in classes else ins or ('DEL-junction' if classes else 'non-SV')
    return dict(mechanism_group='non-SV' if result == 'non-SV' else 'SV-related',
                mechanism_class=result, mechanism_status='evaluated_supplied_ins_del')
