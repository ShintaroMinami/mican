#include "mican.h"
#include "fsscanf.h"
#include <ctype.h>
#include <strings.h>

/*
 * Minimal mmCIF reader.
 *
 * It reads the first _atom_site loop of the file, keeps the rows of the first
 * model only (like stopping at the first ENDMDL in a PDB file), and converts
 * each row to an ATOM_REC that follows PDB field conventions, so that the
 * rest of the program cannot tell the formats apart:
 *   - author items (auth_*) are used when present, label_* otherwise;
 *   - atom names are aligned as in PDB columns 13-16 using the element
 *     symbol (" CA " for C-alpha, "CA  " for calcium);
 *   - residue names are right-justified to 3 characters ("  A" for RNA A);
 *   - "?" and "." (unknown / not applicable) read as a blank PDB field.
 * Supported syntax: data_ blocks, loop_, tags, unquoted, 'single' and
 * "double" quoted values, ;-delimited text fields, and # comments.
 */

#define FIELD_LEN 32

/* token types */
enum { TOK_EOF, TOK_VALUE, TOK_TAG, TOK_LOOP, TOK_DATA, TOK_KEYWORD };

typedef struct {
  FILE *fp;
  char *line;      /* current line, newline stripped */
  size_t cap;
  char *pos;       /* next character to scan in line (NULL: read a line) */
  char *text;      /* contents of the last ;-delimited text field */
  size_t tcap;
  const char *tok; /* last token, valid until the next call */
} LEXER;

/* _atom_site items used */
enum {
  F_GROUP, F_ID, F_TYPE, F_LATOM, F_AATOM, F_LCOMP, F_ACOMP, F_LASYM, F_AASYM,
  F_LSEQ, F_ASEQ, F_X, F_Y, F_Z, F_OCC, F_BISO, F_MODEL, N_FIELD
};
PRIVATE const char *field_name[N_FIELD] = {
    "group_PDB",     "id",           "type_symbol",    "label_atom_id",
    "auth_atom_id",  "label_comp_id", "auth_comp_id",  "label_asym_id",
    "auth_asym_id",  "label_seq_id", "auth_seq_id",    "Cartn_x",
    "Cartn_y",       "Cartn_z",      "occupancy",      "B_iso_or_equiv",
    "pdbx_PDB_model_num"};

PRIVATE void *xrealloc(void *p, size_t size) {
  p = realloc(p, size);
  if (p == NULL) {
    fprintf(stderr, "!! Error   : memory allocation failed\n");
    exit(1);
  }
  return (p);
}

/* Reads a whole line of any length into *buf; strips "\n" / "\r\n". */
PRIVATE int read_line(FILE *fp, char **buf, size_t *cap) {
  size_t n = 0;
  if (*cap == 0) {
    *cap = 512;
    *buf = (char *)xrealloc(NULL, *cap);
  }
  while (fgets(*buf + n, (int)(*cap - n), fp) != NULL) {
    n += strlen(*buf + n);
    if (n == 0 || (*buf)[n - 1] == '\n' || n + 1 < *cap) {
      break; /* end of line, or end of file without a newline */
    }
    *cap *= 2;
    *buf = (char *)xrealloc(*buf, *cap);
  }
  if (n == 0) {
    return (FALSE);
  }
  while (n > 0 && ((*buf)[n - 1] == '\n' || (*buf)[n - 1] == '\r')) {
    (*buf)[--n] = '\0';
  }
  return (TRUE);
}

/* Reads a ;-delimited text field starting at the current line. */
PRIVATE int read_text_field(LEXER *lx) {
  size_t n = 0, len;
  const char *src = lx->line + 1;
  for (;;) {
    len = strlen(src);
    if (n + len + 2 > lx->tcap) {
      lx->tcap = (n + len + 2) * 2;
      lx->text = (char *)xrealloc(lx->text, lx->tcap);
    }
    if (n > 0) {
      lx->text[n++] = '\n';
    }
    memcpy(lx->text + n, src, len + 1);
    n += len;
    if (!read_line(lx->fp, &lx->line, &lx->cap)) {
      lx->pos = NULL; /* unterminated: the field runs to the end of file */
      break;
    }
    if (lx->line[0] == ';') {
      lx->pos = lx->line + 1;
      break;
    }
    src = lx->line;
  }
  lx->tok = lx->text;
  return (TOK_VALUE);
}

PRIVATE int next_token(LEXER *lx) {
  char *p, *start, q;
  for (;;) {
    if (lx->pos == NULL) {
      if (!read_line(lx->fp, &lx->line, &lx->cap)) {
        return (TOK_EOF);
      }
      if (lx->line[0] == ';') {
        return (read_text_field(lx));
      }
      lx->pos = lx->line;
    }
    p = lx->pos;
    while (*p != '\0' && isspace((unsigned char)*p)) {
      p++;
    }
    if (*p == '\0' || *p == '#') {
      lx->pos = NULL;
      continue;
    }
    if (*p == '\'' || *p == '"') {
      /* a quote closes only when followed by whitespace or the line end */
      q = *p;
      start = ++p;
      while (*p != '\0' &&
             !(*p == q && (p[1] == '\0' || isspace((unsigned char)p[1])))) {
        p++;
      }
      lx->pos = (*p != '\0') ? p + 1 : p;
      *p = '\0';
      lx->tok = start;
      return (TOK_VALUE);
    }
    start = p;
    while (*p != '\0' && !isspace((unsigned char)*p)) {
      p++;
    }
    lx->pos = (*p != '\0') ? p + 1 : p;
    *p = '\0';
    lx->tok = start;
    if (start[0] == '_') {
      return (TOK_TAG);
    }
    if (strcasecmp(start, "loop_") == 0) {
      return (TOK_LOOP);
    }
    if (strncasecmp(start, "data_", 5) == 0) {
      return (TOK_DATA);
    }
    if (strncasecmp(start, "save_", 5) == 0 ||
        strcasecmp(start, "global_") == 0 || strcasecmp(start, "stop_") == 0) {
      return (TOK_KEYWORD);
    }
    return (TOK_VALUE);
  }
}

PRIVATE int is_null(const char *v) {
  return (v[0] == '\0' || ((v[0] == '?' || v[0] == '.') && v[1] == '\0'));
}

/* Value of item f in the current row, or NULL if absent or unknown. */
PRIVATE const char *item(char (*val)[FIELD_LEN], const int *have, int f) {
  return ((have[f] && !is_null(val[f])) ? val[f] : NULL);
}

/* Author item if given, label item otherwise. */
PRIVATE const char *item2(char (*val)[FIELD_LEN], const int *have, int f_auth,
                          int f_label) {
  const char *v = item(val, have, f_auth);
  return ((v != NULL) ? v : item(val, have, f_label));
}

PRIVATE int to_int(const char *v) {
  int x = 0;
  if (v != NULL) {
    scan_int(v, &x);
  }
  return (x);
}

PRIVATE float to_float(const char *v) {
  float x = (float)0.0;
  if (v != NULL) {
    scan_float(v, &x);
  }
  return (x);
}

PRIVATE void row_to_atom(char (*val)[FIELD_LEN], const int *have,
                         ATOM_REC *a) {
  const char *group = item(val, have, F_GROUP);
  const char *atom = item2(val, have, F_AATOM, F_LATOM);
  const char *element = item(val, have, F_TYPE);
  const char *comp = item2(val, have, F_ACOMP, F_LCOMP);
  const char *asym = item2(val, have, F_AASYM, F_LASYM);

  snprintf(a->record, sizeof(a->record), "%-6.6s",
           (group != NULL) ? group : "ATOM");
  a->serial = to_int(item(val, have, F_ID));
  if (atom == NULL) {
    atom = "";
  }
  /* PDB alignment: a name shorter than 4 characters starts in column 14,
     unless its element symbol has two letters */
  if (strlen(atom) >= 4 || (element != NULL && strlen(element) == 2)) {
    snprintf(a->name, sizeof(a->name), "%-4.4s", atom);
  } else {
    snprintf(a->name, sizeof(a->name), " %-3s", atom);
  }
  snprintf(a->resName, sizeof(a->resName), "%3.3s",
           (comp != NULL) ? comp : "");
  snprintf(a->chain, sizeof(a->chain), "%.4s", (asym != NULL) ? asym : " ");
  a->resSeq = to_int(item2(val, have, F_ASEQ, F_LSEQ));
  a->x = to_float(item(val, have, F_X));
  a->y = to_float(item(val, have, F_Y));
  a->z = to_float(item(val, have, F_Z));
  a->occupancy = to_float(item(val, have, F_OCC));
  a->tempFactor = to_float(item(val, have, F_BISO));
}

/*************************************/
/*    FUNCTION  is_mmcif_file        */
/*************************************/
/* A file is mmCIF when its first line that is neither blank nor a comment
   starts with "data_". */
int is_mmcif_file(const char *filename) {
  FILE *fp;
  char *line = NULL, *p;
  size_t cap = 0;
  int cif = FALSE;

  if ((fp = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "!! Error   : No such file (%s)\n", filename);
    exit(1);
  }
  while (read_line(fp, &line, &cap)) {
    p = line;
    while (*p != '\0' && isspace((unsigned char)*p)) {
      p++;
    }
    if (*p == '\0' || *p == '#') {
      continue;
    }
    cif = (strncasecmp(p, "data_", 5) == 0);
    break;
  }
  free(line);
  fclose(fp);
  return (cif);
}

/*************************************/
/*    FUNCTION  read_mmcif_atoms     */
/*************************************/
/* Reads the atoms of the first model into a new array *atoms (0-based) and
   returns their number. */
int read_mmcif_atoms(const char *filename, ATOM_REC **atoms) {
  LEXER lx = {NULL, NULL, 0, NULL, NULL, 0, NULL};
  int tok, ncol = 0, ccap = 0, icol, f, natom = 0, acap = 0;
  int *colmap = NULL; /* loop column -> item index (-1: not used) */
  int have[N_FIELD];
  char val[N_FIELD][FIELD_LEN];
  char model0[FIELD_LEN] = "";
  int found = FALSE, stop = FALSE;
  ATOM_REC *a = NULL;

  if ((lx.fp = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "!! Error   : No such file (%s)\n", filename);
    exit(1);
  }

  tok = next_token(&lx);
  while (tok != TOK_EOF && !found) {
    if (tok != TOK_LOOP) {
      tok = next_token(&lx);
      continue;
    }
    /** loop header: the first tag decides the category **/
    ncol = 0;
    while ((tok = next_token(&lx)) == TOK_TAG) {
      if (ncol == 0) {
        found = (strncasecmp(lx.tok, "_atom_site.", 11) == 0);
      }
      if (ncol >= ccap) {
        ccap = 2 * ccap + 16;
        colmap = (int *)xrealloc(colmap, (size_t)ccap * sizeof(int));
      }
      colmap[ncol] = -1;
      for (f = 0; found && strncasecmp(lx.tok, "_atom_site.", 11) == 0 &&
                  f < N_FIELD;
           f++) {
        if (strcasecmp(lx.tok + 11, field_name[f]) == 0) {
          colmap[ncol] = f;
        }
      }
      ncol++;
    }
  }
  if (!found) {
    fprintf(stderr, "!! Error   : No _atom_site loop in mmCIF file (%s)\n",
            filename);
    exit(1);
  }

  /** required items **/
  for (f = 0; f < N_FIELD; f++) {
    have[f] = FALSE;
    val[f][0] = '\0';
  }
  for (icol = 0; icol < ncol; icol++) {
    if (colmap[icol] >= 0) {
      have[colmap[icol]] = TRUE;
    }
  }
  if (!have[F_X] || !have[F_Y] || !have[F_Z] ||
      !(have[F_LATOM] || have[F_AATOM]) || !(have[F_LCOMP] || have[F_ACOMP]) ||
      !(have[F_LASYM] || have[F_AASYM]) || !(have[F_LSEQ] || have[F_ASEQ])) {
    fprintf(stderr,
            "!! Error   : _atom_site needs Cartn_x/y/z and atom_id, comp_id, "
            "asym_id, seq_id (label or auth) (%s)\n",
            filename);
    exit(1);
  }

  /** rows **/
  icol = 0;
  for (; tok == TOK_VALUE && !stop; tok = next_token(&lx)) {
    if (colmap[icol] >= 0) {
      snprintf(val[colmap[icol]], FIELD_LEN, "%s", lx.tok);
    }
    if (++icol < ncol) {
      continue;
    }
    icol = 0;
    /* first model only */
    if (have[F_MODEL]) {
      if (natom == 0) {
        strcpy(model0, val[F_MODEL]);
      } else if (strcmp(val[F_MODEL], model0) != 0) {
        stop = TRUE;
        continue;
      }
    }
    if (natom >= acap) {
      acap = 2 * acap + 1024;
      a = (ATOM_REC *)xrealloc(a, (size_t)acap * sizeof(ATOM_REC));
    }
    row_to_atom(val, have, &a[natom]);
    natom++;
  }

  free(colmap);
  free(lx.line);
  free(lx.text);
  fclose(lx.fp);
  *atoms = a;
  return (natom);
}
