# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: ga4gh/rna_quantification_service.proto

import sys
_b=sys.version_info[0]<3 and (lambda x:x) or (lambda x:x.encode('latin1'))
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
from google.protobuf import descriptor_pb2
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


from ga4gh import rna_quantification_pb2 as ga4gh_dot_rna__quantification__pb2


DESCRIPTOR = _descriptor.FileDescriptor(
  name='ga4gh/rna_quantification_service.proto',
  package='ga4gh',
  syntax='proto3',
  serialized_pb=_b('\n&ga4gh/rna_quantification_service.proto\x12\x05ga4gh\x1a\x1ega4gh/rna_quantification.proto\"_\n\"SearchRnaQuantificationSetsRequest\x12\x12\n\ndataset_id\x18\x01 \x01(\t\x12\x11\n\tpage_size\x18\x02 \x01(\x05\x12\x12\n\npage_token\x18\x03 \x01(\t\"|\n#SearchRnaQuantificationSetsResponse\x12<\n\x17rna_quantification_sets\x18\x01 \x03(\x0b\x32\x1b.ga4gh.RnaQuantificationSet\x12\x17\n\x0fnext_page_token\x18\x02 \x01(\t\"C\n\x1eGetRnaQuantificationSetRequest\x12!\n\x19rna_quantification_set_id\x18\x01 \x01(\t\"\x7f\n\x1fSearchRnaQuantificationsRequest\x12!\n\x19rna_quantification_set_id\x18\x01 \x01(\t\x12\x12\n\ndataset_id\x18\x02 \x01(\t\x12\x11\n\tpage_size\x18\x03 \x01(\x05\x12\x12\n\npage_token\x18\x04 \x01(\t\"r\n SearchRnaQuantificationsResponse\x12\x35\n\x13rna_quantifications\x18\x01 \x03(\x0b\x32\x18.ga4gh.RnaQuantification\x12\x17\n\x0fnext_page_token\x18\x02 \x01(\t\"<\n\x1bGetRnaQuantificationRequest\x12\x1d\n\x15rna_quantification_id\x18\x01 \x01(\t\"\x92\x01\n\x1dSearchExpressionLevelsRequest\x12\x18\n\x10\x66\x65\x61ture_group_id\x18\x01 \x01(\t\x12\x1d\n\x15rna_quantification_id\x18\x02 \x01(\t\x12\x11\n\tthreshold\x18\x03 \x01(\x02\x12\x11\n\tpage_size\x18\x04 \x01(\x05\x12\x12\n\npage_token\x18\x05 \x01(\t\"l\n\x1eSearchExpressionLevelsResponse\x12\x31\n\x11\x65xpression_levels\x18\x01 \x03(\x0b\x32\x16.ga4gh.ExpressionLevel\x12\x17\n\x0fnext_page_token\x18\x02 \x01(\t\"8\n\x19GetExpressionLevelRequest\x12\x1b\n\x13\x65xpression_level_id\x18\x01 \x01(\t\"W\n\x1aSearchFeatureGroupsRequest\x12\x12\n\ndataset_id\x18\x01 \x01(\t\x12\x11\n\tpage_size\x18\x02 \x01(\x05\x12\x12\n\npage_token\x18\x03 \x01(\t\"c\n\x1bSearchFeatureGroupsResponse\x12+\n\x0e\x66\x65\x61ture_groups\x18\x01 \x03(\x0b\x32\x13.ga4gh.FeatureGroup\x12\x17\n\x0fnext_page_token\x18\x02 \x01(\t\"2\n\x16GetFeatureGroupRequest\x12\x18\n\x10\x66\x65\x61ture_group_id\x18\x01 \x01(\t2\x8e\x06\n\x18RnaQuantificationService\x12t\n\x1bSearchRnaQuantificationSets\x12).ga4gh.SearchRnaQuantificationSetsRequest\x1a*.ga4gh.SearchRnaQuantificationSetsResponse\x12]\n\x17GetRnaQuantificationSet\x12%.ga4gh.GetRnaQuantificationSetRequest\x1a\x1b.ga4gh.RnaQuantificationSet\x12k\n\x18SearchRnaQuantifications\x12&.ga4gh.SearchRnaQuantificationsRequest\x1a\'.ga4gh.SearchRnaQuantificationsResponse\x12T\n\x14GetRnaQuantification\x12\".ga4gh.GetRnaQuantificationRequest\x1a\x18.ga4gh.RnaQuantification\x12\x65\n\x16SearchExpressionLevels\x12$.ga4gh.SearchExpressionLevelsRequest\x1a%.ga4gh.SearchExpressionLevelsResponse\x12N\n\x12GetExpressionLevel\x12 .ga4gh.GetExpressionLevelRequest\x1a\x16.ga4gh.ExpressionLevel\x12\\\n\x13SearchFeatureGroups\x12!.ga4gh.SearchFeatureGroupsRequest\x1a\".ga4gh.SearchFeatureGroupsResponse\x12\x45\n\x0fGetFeatureGroup\x12\x1d.ga4gh.GetFeatureGroupRequest\x1a\x13.ga4gh.FeatureGroupb\x06proto3')
  ,
  dependencies=[ga4gh_dot_rna__quantification__pb2.DESCRIPTOR,])
_sym_db.RegisterFileDescriptor(DESCRIPTOR)




_SEARCHRNAQUANTIFICATIONSETSREQUEST = _descriptor.Descriptor(
  name='SearchRnaQuantificationSetsRequest',
  full_name='ga4gh.SearchRnaQuantificationSetsRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='dataset_id', full_name='ga4gh.SearchRnaQuantificationSetsRequest.dataset_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_size', full_name='ga4gh.SearchRnaQuantificationSetsRequest.page_size', index=1,
      number=2, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_token', full_name='ga4gh.SearchRnaQuantificationSetsRequest.page_token', index=2,
      number=3, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=81,
  serialized_end=176,
)


_SEARCHRNAQUANTIFICATIONSETSRESPONSE = _descriptor.Descriptor(
  name='SearchRnaQuantificationSetsResponse',
  full_name='ga4gh.SearchRnaQuantificationSetsResponse',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='rna_quantification_sets', full_name='ga4gh.SearchRnaQuantificationSetsResponse.rna_quantification_sets', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='next_page_token', full_name='ga4gh.SearchRnaQuantificationSetsResponse.next_page_token', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=178,
  serialized_end=302,
)


_GETRNAQUANTIFICATIONSETREQUEST = _descriptor.Descriptor(
  name='GetRnaQuantificationSetRequest',
  full_name='ga4gh.GetRnaQuantificationSetRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='rna_quantification_set_id', full_name='ga4gh.GetRnaQuantificationSetRequest.rna_quantification_set_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=304,
  serialized_end=371,
)


_SEARCHRNAQUANTIFICATIONSREQUEST = _descriptor.Descriptor(
  name='SearchRnaQuantificationsRequest',
  full_name='ga4gh.SearchRnaQuantificationsRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='rna_quantification_set_id', full_name='ga4gh.SearchRnaQuantificationsRequest.rna_quantification_set_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='dataset_id', full_name='ga4gh.SearchRnaQuantificationsRequest.dataset_id', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_size', full_name='ga4gh.SearchRnaQuantificationsRequest.page_size', index=2,
      number=3, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_token', full_name='ga4gh.SearchRnaQuantificationsRequest.page_token', index=3,
      number=4, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=373,
  serialized_end=500,
)


_SEARCHRNAQUANTIFICATIONSRESPONSE = _descriptor.Descriptor(
  name='SearchRnaQuantificationsResponse',
  full_name='ga4gh.SearchRnaQuantificationsResponse',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='rna_quantifications', full_name='ga4gh.SearchRnaQuantificationsResponse.rna_quantifications', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='next_page_token', full_name='ga4gh.SearchRnaQuantificationsResponse.next_page_token', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=502,
  serialized_end=616,
)


_GETRNAQUANTIFICATIONREQUEST = _descriptor.Descriptor(
  name='GetRnaQuantificationRequest',
  full_name='ga4gh.GetRnaQuantificationRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='rna_quantification_id', full_name='ga4gh.GetRnaQuantificationRequest.rna_quantification_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=618,
  serialized_end=678,
)


_SEARCHEXPRESSIONLEVELSREQUEST = _descriptor.Descriptor(
  name='SearchExpressionLevelsRequest',
  full_name='ga4gh.SearchExpressionLevelsRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='feature_group_id', full_name='ga4gh.SearchExpressionLevelsRequest.feature_group_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='rna_quantification_id', full_name='ga4gh.SearchExpressionLevelsRequest.rna_quantification_id', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='threshold', full_name='ga4gh.SearchExpressionLevelsRequest.threshold', index=2,
      number=3, type=2, cpp_type=6, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_size', full_name='ga4gh.SearchExpressionLevelsRequest.page_size', index=3,
      number=4, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_token', full_name='ga4gh.SearchExpressionLevelsRequest.page_token', index=4,
      number=5, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=681,
  serialized_end=827,
)


_SEARCHEXPRESSIONLEVELSRESPONSE = _descriptor.Descriptor(
  name='SearchExpressionLevelsResponse',
  full_name='ga4gh.SearchExpressionLevelsResponse',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='expression_levels', full_name='ga4gh.SearchExpressionLevelsResponse.expression_levels', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='next_page_token', full_name='ga4gh.SearchExpressionLevelsResponse.next_page_token', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=829,
  serialized_end=937,
)


_GETEXPRESSIONLEVELREQUEST = _descriptor.Descriptor(
  name='GetExpressionLevelRequest',
  full_name='ga4gh.GetExpressionLevelRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='expression_level_id', full_name='ga4gh.GetExpressionLevelRequest.expression_level_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=939,
  serialized_end=995,
)


_SEARCHFEATUREGROUPSREQUEST = _descriptor.Descriptor(
  name='SearchFeatureGroupsRequest',
  full_name='ga4gh.SearchFeatureGroupsRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='dataset_id', full_name='ga4gh.SearchFeatureGroupsRequest.dataset_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_size', full_name='ga4gh.SearchFeatureGroupsRequest.page_size', index=1,
      number=2, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_token', full_name='ga4gh.SearchFeatureGroupsRequest.page_token', index=2,
      number=3, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=997,
  serialized_end=1084,
)


_SEARCHFEATUREGROUPSRESPONSE = _descriptor.Descriptor(
  name='SearchFeatureGroupsResponse',
  full_name='ga4gh.SearchFeatureGroupsResponse',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='feature_groups', full_name='ga4gh.SearchFeatureGroupsResponse.feature_groups', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='next_page_token', full_name='ga4gh.SearchFeatureGroupsResponse.next_page_token', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=1086,
  serialized_end=1185,
)


_GETFEATUREGROUPREQUEST = _descriptor.Descriptor(
  name='GetFeatureGroupRequest',
  full_name='ga4gh.GetFeatureGroupRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='feature_group_id', full_name='ga4gh.GetFeatureGroupRequest.feature_group_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=1187,
  serialized_end=1237,
)

_SEARCHRNAQUANTIFICATIONSETSRESPONSE.fields_by_name['rna_quantification_sets'].message_type = ga4gh_dot_rna__quantification__pb2._RNAQUANTIFICATIONSET
_SEARCHRNAQUANTIFICATIONSRESPONSE.fields_by_name['rna_quantifications'].message_type = ga4gh_dot_rna__quantification__pb2._RNAQUANTIFICATION
_SEARCHEXPRESSIONLEVELSRESPONSE.fields_by_name['expression_levels'].message_type = ga4gh_dot_rna__quantification__pb2._EXPRESSIONLEVEL
_SEARCHFEATUREGROUPSRESPONSE.fields_by_name['feature_groups'].message_type = ga4gh_dot_rna__quantification__pb2._FEATUREGROUP
DESCRIPTOR.message_types_by_name['SearchRnaQuantificationSetsRequest'] = _SEARCHRNAQUANTIFICATIONSETSREQUEST
DESCRIPTOR.message_types_by_name['SearchRnaQuantificationSetsResponse'] = _SEARCHRNAQUANTIFICATIONSETSRESPONSE
DESCRIPTOR.message_types_by_name['GetRnaQuantificationSetRequest'] = _GETRNAQUANTIFICATIONSETREQUEST
DESCRIPTOR.message_types_by_name['SearchRnaQuantificationsRequest'] = _SEARCHRNAQUANTIFICATIONSREQUEST
DESCRIPTOR.message_types_by_name['SearchRnaQuantificationsResponse'] = _SEARCHRNAQUANTIFICATIONSRESPONSE
DESCRIPTOR.message_types_by_name['GetRnaQuantificationRequest'] = _GETRNAQUANTIFICATIONREQUEST
DESCRIPTOR.message_types_by_name['SearchExpressionLevelsRequest'] = _SEARCHEXPRESSIONLEVELSREQUEST
DESCRIPTOR.message_types_by_name['SearchExpressionLevelsResponse'] = _SEARCHEXPRESSIONLEVELSRESPONSE
DESCRIPTOR.message_types_by_name['GetExpressionLevelRequest'] = _GETEXPRESSIONLEVELREQUEST
DESCRIPTOR.message_types_by_name['SearchFeatureGroupsRequest'] = _SEARCHFEATUREGROUPSREQUEST
DESCRIPTOR.message_types_by_name['SearchFeatureGroupsResponse'] = _SEARCHFEATUREGROUPSRESPONSE
DESCRIPTOR.message_types_by_name['GetFeatureGroupRequest'] = _GETFEATUREGROUPREQUEST

SearchRnaQuantificationSetsRequest = _reflection.GeneratedProtocolMessageType('SearchRnaQuantificationSetsRequest', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHRNAQUANTIFICATIONSETSREQUEST,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchRnaQuantificationSetsRequest)
  ))
_sym_db.RegisterMessage(SearchRnaQuantificationSetsRequest)

SearchRnaQuantificationSetsResponse = _reflection.GeneratedProtocolMessageType('SearchRnaQuantificationSetsResponse', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHRNAQUANTIFICATIONSETSRESPONSE,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchRnaQuantificationSetsResponse)
  ))
_sym_db.RegisterMessage(SearchRnaQuantificationSetsResponse)

GetRnaQuantificationSetRequest = _reflection.GeneratedProtocolMessageType('GetRnaQuantificationSetRequest', (_message.Message,), dict(
  DESCRIPTOR = _GETRNAQUANTIFICATIONSETREQUEST,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.GetRnaQuantificationSetRequest)
  ))
_sym_db.RegisterMessage(GetRnaQuantificationSetRequest)

SearchRnaQuantificationsRequest = _reflection.GeneratedProtocolMessageType('SearchRnaQuantificationsRequest', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHRNAQUANTIFICATIONSREQUEST,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchRnaQuantificationsRequest)
  ))
_sym_db.RegisterMessage(SearchRnaQuantificationsRequest)

SearchRnaQuantificationsResponse = _reflection.GeneratedProtocolMessageType('SearchRnaQuantificationsResponse', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHRNAQUANTIFICATIONSRESPONSE,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchRnaQuantificationsResponse)
  ))
_sym_db.RegisterMessage(SearchRnaQuantificationsResponse)

GetRnaQuantificationRequest = _reflection.GeneratedProtocolMessageType('GetRnaQuantificationRequest', (_message.Message,), dict(
  DESCRIPTOR = _GETRNAQUANTIFICATIONREQUEST,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.GetRnaQuantificationRequest)
  ))
_sym_db.RegisterMessage(GetRnaQuantificationRequest)

SearchExpressionLevelsRequest = _reflection.GeneratedProtocolMessageType('SearchExpressionLevelsRequest', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHEXPRESSIONLEVELSREQUEST,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchExpressionLevelsRequest)
  ))
_sym_db.RegisterMessage(SearchExpressionLevelsRequest)

SearchExpressionLevelsResponse = _reflection.GeneratedProtocolMessageType('SearchExpressionLevelsResponse', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHEXPRESSIONLEVELSRESPONSE,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchExpressionLevelsResponse)
  ))
_sym_db.RegisterMessage(SearchExpressionLevelsResponse)

GetExpressionLevelRequest = _reflection.GeneratedProtocolMessageType('GetExpressionLevelRequest', (_message.Message,), dict(
  DESCRIPTOR = _GETEXPRESSIONLEVELREQUEST,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.GetExpressionLevelRequest)
  ))
_sym_db.RegisterMessage(GetExpressionLevelRequest)

SearchFeatureGroupsRequest = _reflection.GeneratedProtocolMessageType('SearchFeatureGroupsRequest', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHFEATUREGROUPSREQUEST,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchFeatureGroupsRequest)
  ))
_sym_db.RegisterMessage(SearchFeatureGroupsRequest)

SearchFeatureGroupsResponse = _reflection.GeneratedProtocolMessageType('SearchFeatureGroupsResponse', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHFEATUREGROUPSRESPONSE,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchFeatureGroupsResponse)
  ))
_sym_db.RegisterMessage(SearchFeatureGroupsResponse)

GetFeatureGroupRequest = _reflection.GeneratedProtocolMessageType('GetFeatureGroupRequest', (_message.Message,), dict(
  DESCRIPTOR = _GETFEATUREGROUPREQUEST,
  __module__ = 'ga4gh.rna_quantification_service_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.GetFeatureGroupRequest)
  ))
_sym_db.RegisterMessage(GetFeatureGroupRequest)


# @@protoc_insertion_point(module_scope)