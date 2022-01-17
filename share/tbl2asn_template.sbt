Submit-block ::= {
  contact {
    contact {
      name name {
        last "[% last_name %]",
        first "[% first_name %]",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "[% affiliation %]",
        div "[% division %]",
        city "[% city %]",
        sub "[% state %]",
        country "[% country %]",
        street "[% street %]",
        email "[% email %]",
        postal-code "[% zip %]"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "[% last_name %]",
            first "[% first_name %]",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "[% affiliation %]",
        div "[% division %]",
        city "[% city %]",
        sub "[% state %]",
        country "[% country %]",
        street "[% street %]",
        postal-code "[% zip %]"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "[% last_name %]",
              first "[% first_name %]",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "[% hit_description %]"
    }
  }
}
Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        "[% bioproject %]"
      }
    },
    {
      label str "BioSample",
      num 1,
      data strs {
        "[% biosample %]"
      }
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "[% taxon %]",
      data str "ALT EMAIL:[% second_email %]"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "[% user_comment %]",
      data str "[% hit_description %]"
    }
  }
}
