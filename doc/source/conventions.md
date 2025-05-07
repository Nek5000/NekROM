# Coding Conventions for MOR

1. Capitalization

Lowercase characters should be used everywhere except for names of included files.

2. Indentation

Each Indentation should be 3 characters wide.

3. White-space

No trailing white-spaces. No blank comment lines.

4. Subroutines

Brevity in subroutine naming is preferred, but not at the expense of clarity. Ideally, subroutines should fit in one screen.

5. Variables

Variable names must be short, but common block variables should be descriptive. Vector variables in general should not be allocated inside a subroutine. Instead, they should be declared in a common block.

6. Comments

Brief description should be provided for non-trivial subroutines. The input/output arguments should be described as well.

7. Common Blocks

Limit use of common blocks i.e., functional subroutines are preferred.

8. For-Loops

Use do-enddo delimiters for for-loops. Directly nested for-loops should not be indented:

```
do k=k0,k1
do j=j0,j1
do i=i0,i1
   ! body here (i,j,k)
enddo
enddo
enddo
```
