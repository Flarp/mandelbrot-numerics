#include <mandelbrot-numerics.h>

extern void m_d_mat2_set(m_d_mat2 *o, const m_d_mat2 *m) {
  o->a = m->a;
  o->b = m->b;
  o->c = m->c;
  o->d = m->d;
}

extern void m_d_mat2_id(m_d_mat2 *o) {
  o->a = 1;
  o->b = 0;
  o->c = 0;
  o->d = 1;
}

extern double _Complex m_d_mat2_tr(const m_d_mat2 *m) {
  return m->a + m->d;
}

extern double _Complex m_d_mat2_det(const m_d_mat2 *m) {
  return m->a * m->d - m->b * m->c;
}

extern void m_d_mat2_inv(m_d_mat2 *m1, const m_d_mat2 *m) {
  double _Complex det = m_d_mat2_det(m);
  double _Complex a =  m->d / det;
  double _Complex b = -m->b / det;
  double _Complex c = -m->c / det;
  double _Complex d =  m->a / det;
  m1->a = a;
  m1->b = b;
  m1->c = c;
  m1->d = d;
}

extern void m_d_mat2_mul(m_d_mat2 *o, const m_d_mat2 *l, const m_d_mat2 *r) {
  double _Complex a, b, c, d;
  a = l->a * r->a + l->b * r->c;
  b = l->a * r->b + l->b * r->d;
  c = l->c * r->a + l->d * r->c;
  d = l->c * r->b + l->d * r->d;
  o->a = a;
  o->b = b;
  o->c = c;
  o->d = d;
}

extern void m_d_mat2_diagonalize(m_d_mat2 *p, m_d_mat2 *d, m_d_mat2 *p1, const m_d_mat2 *m) {
  double _Complex tr2 = m_d_mat2_tr(m) / 2;
  double _Complex det = m_d_mat2_det(m);
  double _Complex k = csqrt(tr2 * tr2 - det);
  double _Complex l1 = tr2 + k;
  double _Complex l2 = tr2 - k;
  d->a = l1;
  d->b = 0;
  d->c = 0;
  d->d = l2;
  if (m->b != 0) {
    p->a = m->b;
    p->b = m->b;
    p->c = l1 - m->a;
    p->d = l2 - m->a;
    m_d_mat2_inv(p1, p);
  } else if (m->c != 0) {
    p->a = l1 - m->d;
    p->b = l2 - m->d;
    p->c = m->c;
    p->d = m->c;
    m_d_mat2_inv(p1, p);
  } else {
    p1->a = p->a = 1;
    p1->b = p->b = 0;
    p1->c = p->c = 0;
    p1->d = p->d = 1;
  }
}

extern void m_d_mat2_moebius3(m_d_mat2 *m, double _Complex zero, double _Complex one, double _Complex infinity) {
  m->a = infinity * (zero - one);
  m->b = zero * (one - infinity);
  m->c = zero - one;
  m->d = one - infinity;
}

struct m_d_mat2_interp {
  m_d_mat2 fp, d, p1;
};

extern m_d_mat2_interp *m_d_mat2_interp_new(void) {
  return malloc(sizeof(m_d_mat2_interp));
}

extern void m_d_mat2_interp_delete(m_d_mat2_interp *i) {
  free(i);
}

extern void m_d_mat2_interp_init(m_d_mat2_interp *i, const m_d_mat2 *f, const m_d_mat2 *g) {
  m_d_mat2 f1;
  m_d_mat2_inv(&f1, f);
  m_d_mat2 f1g;
  m_d_mat2_mul(&f1g, &f1, g);
  m_d_mat2 p;
  m_d_mat2_diagonalize(&p, &i->d, &i->p1, &f1g);
  m_d_mat2_mul(&i->fp, f, &p);
}

extern void m_d_mat2_interp_do(m_d_mat2 *m, const m_d_mat2_interp *i, double t) {
  m_d_mat2 e;
  e.a = cpow(i->d.a, t);
  e.b = 0;
  e.c = 0;
  e.d = cpow(i->d.d, t);
  m_d_mat2 fpe;
  m_d_mat2_mul(&fpe, &i->fp, &e);
  m_d_mat2_mul(m, &fpe, &i->p1);
}
